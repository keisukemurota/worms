import sys
import numpy as np
import os
import logging
import re
import pandas as pd
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import utils  # noqa: E402
from utils.parser import get_parser  # noqa: E402
from utils import get_logger, extract_info_from_file  # noqa: E402


parser = get_parser()
logger = get_logger("log.log", stdout=True, level=logging.INFO)

parser.add_argument(
    '-f',
    '--path',
    type=str,
    required=True,
    help='The python file will look for the Hamiltonian and the best unitary in this directory.')

parser.add_argument('-s', '--sweeps', type=int, required=True,
                    help='The number of sweeps to perform.')
parser.add_argument('--original', action='store_true', default=False,
                    help='If this flag is set, the original Hamiltonian will be used.')

parser.add_argument('-k', '--top_k', type=int, default=1,
                    help='The number of top unitary paths to consider.')
args = parser.parse_args()
search_path = Path(args.path)

# check if given path exists
if not os.path.isdir(args.path):
    logging.debug("current dir is: ", os.getcwd())
    logging.info("given serch path doesn't exit.")
    logging.info("search path: {}".format(args.path))
    sys.exit()

rmsKit_directory = Path(os.path.abspath(__file__)).parent.parent


if __name__ == "__main__":

    # define parameters list to be passed to the run_worm function
    # beta = np.linspace(0.5, 5, 3)
    if args.model == "SS2D":
        beta = np.array([0.25, 1, 4, 10])
        # beta = np.array([10])
        L_list = [[2, 2], [3, 3]]
        logger.info("RUN SS2D MODEL")
    elif args.model == "HXYZ2D":
        beta = np.array([0.5, 1, 4])
        L_list = [[3, 3], [4, 4]]
        logger.info("RUN HXYZ2D MODEL")
    elif args.model == "KH2D":
        beta = np.array([1, 4])
        L_list = [[4, 4], [5, 5]]
        logger.info("RUN HXYZ2D MODEL")
    elif args.model == "BLBQ1D":
        beta = np.array([1, 4])
        L_list = [[10], [11]]
        beta_select = 1
        L_list_select = [4]
        logger.info("RUN BLBQ1D MODEL")
    elif args.model == "FF2D":
        beta = np.array([1])
        L_list = [[4, 4]]
        logger.info("RUN FF2D MODEL")
    else:
        beta = np.array([1, 4])
        L_list = [[10], [11]]
        logger.info("RUN {} MODEL".format(args.model))

    T_list = 1/beta
    # define the number of samples
    p = args.num_threads
    M = args.sweeps

    if M % p != 0:
        logger.info(
            "Warning: M is not divisible by p. The number of sweeps will be rounded down.")
        M = (M // p) * p

    # min_path, min_loss, ham_path = utils.path_with_lowest_loss(
    #     args.path, return_ham=True, absolute_path=True)
    search_path = Path(args.path)
    if search_path.is_symlink():
        search_path = search_path.resolve()
        logging.warning("The given path is a symbolic link.")
        logging.warning("The path will be resolved to {}".format(search_path))

    min_loss, init_loss, min_path, ham_path, info_txt_path = utils.path.get_worm_path(
        search_path)

    simu_setting = "sweeps_{}_p_{}".format(M, p)
    if args.original:
        simu_setting += "_original"

    loss_func = info_txt_path.parent.name
    save_path = info_txt_path.parent / "summary" / simu_setting

    # n: create the directory if it does not exist
    save_path.mkdir(parents=True, exist_ok=True)

    # n : find first save_path_{i} that does not exist
    i = 0
    while (save_path / "{}.csv".format(i)).exists():
        i += 1
    save_path = save_path / "{}.csv".format(i)

    logger.info("min_path: {}".format(min_path))
    logger.info("min_loss: {}".format(min_loss))
    logger.info("ham_path: {}".format(ham_path))

    logger.info("model_name: {}".format(args.model))
    logger.info(
        "Search path : {}".format(
            search_path.resolve().as_posix()))
    logger.info("L_list: {}".format(L_list))
    logger.info("T_list: {}".format(T_list))
    logger.info("Summary will be saved to : {}".format(save_path))

    if args.top_k > 1 and beta_select is not None and not args.original:
        logger.info("Getting the top {} unitary paths".format(args.top_k))
        top_k_unitary_paths = utils.path.top_k_upath(search_path, args.top_k)
        logger.info("top_k_unitary_paths: {}".format(top_k_unitary_paths))

        neg_vals = []
        for loss, path in top_k_unitary_paths:
            logger.info("test simulation with loss: {}, path: {}".format(loss, path))
            N_select = 10**5
            N_select_per_cpu = N_select // p
            subprocess_out = utils.run_worm(
                args.model,
                ham_path,
                path,
                L_list_select,
                1/beta_select,
                N_select_per_cpu,
                n=p,
                logging=True,
                obc=args.obc,
                project_dir=rmsKit_directory.parent.parent.resolve())
            output = subprocess_out.stdout.decode("utf-8")

            match = re.search(r'The result will be written in : "(.+?\.txt)"', output)
            try:
                result_file_path = match.group(1)
                data = extract_info_from_file(
                    result_file_path, warning=True, allow_missing=False)
                
                neg_val = data["as_error"] / data["as"]
                neg_vals.append((neg_val, loss, path))
            except Exception as e:
                logger.error("Exception: {}".format(e))
                logger.error(
                    "No result file found. This may be due to an error in the logging in main_MPI.")
                logger.error("subprocess_out: {}".format(output))
                continue
        
        neg_vals.sort(key=lambda x: x[0])

        min_path = neg_vals[0][2]

        logger.info("selected min_path: {} with negativity {}".format(min_path, neg_vals[0][0]))
        logger.info("other negative values: {}".format([(neg, loss) for neg, loss, path in neg_vals[1:]]))



    # run the simulation
    data_list = []
    for L in L_list:
        for T in T_list:
            subprocess_out = utils.run_worm(
                args.model,
                ham_path,
                min_path if not args.original else "",
                L,
                T,
                M,
                n=p,
                logging=True,
                obc=args.obc,
                project_dir=rmsKit_directory.parent.parent.resolve())
            output = subprocess_out.stdout.decode("utf-8")
            # Extract the path using regex
            match = re.search(r'The result will be written in : "(.+?\.txt)"', output)
            try:
                result_file_path = match.group(1)
                logger.info("result_file_path: {}".format(result_file_path))
                data = extract_info_from_file(
                    result_file_path, warning=True, allow_missing=False)

                # n: Store the result in the simu_res.txt file
                data["loss_func"] = loss_func
                data_list.append(data)

                # n: check if the simulation is reliable
                if (data["as_error"] / data["as"]) > 0.2 or (data["as"] <= 0):
                    logger.info(
                        """Negativity was too high {} Simulation is not reliable.
                        The simulation for the following temperature can be ignored.
                        """.format(data["as_error"] / data["as"]))
                else:
                    logger.info(
                        "Simulation succeeded. Sweeps : {} L : {}, T : {}, Negativity : {}".format(
                            M, L, T, data["as_error"] / data["as"]))
            except Exception as e:
                logger.error("Exception: {}".format(e))
                logger.error(
                    "No result file found. This may be due to an error in the logging in main_MPI.")
                logger.error("subprocess_out: {}".format(output))
                continue

    # save the result
    df = pd.DataFrame(data_list)
    df.to_csv(save_path, index=False)
    logging.info("Summary saved to : {}".format(save_path))
