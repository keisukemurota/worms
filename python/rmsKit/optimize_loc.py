"""Optimize the local unitary matrix for the given model to reduce the negative sign problem.

Currently availabe models are:
    1D: HXYZ, BLBQ, Frustration-Free(FF)
    2D: HXYZ, FF

Available loss functions are:
    mel: minimum energy loss. The basic loss function.
    qsmel: quasi-system minimum energy loss. This consider the system hamiltonian with given system size.
    none: no loss function. Just return the hamiltonian
"""
import torch
import numpy as np
import rms_torch
import time
import platform
import logging
from pathlib import Path
import os
import re

from lattice import save_npy, get_model
from utils.parser import get_parser, get_params_parser, get_dtype, print_args_info
from utils import NOW, get_logger, stoquastic

parser = get_parser()

# add symbolic link to generated hamiltonian and unitary
parser.add_argument(
    "--symoblic_link",
    type=Path,  # Changed to Path for direct pathlib support
    default=None,
    help="The symbolic link to the directory the generated hamiltonian and unitary are stored.",
)
args, params, hash_str = get_params_parser(parser)

CHECK = True

if __name__ == "__main__":

    sps = args.sps
    epochs = args.epochs

    log_dir = Path("optimizer_output")
    log_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists
    log_filename = log_dir / f"{NOW}_{hash_str}.log"
    logger = get_logger(log_filename.resolve().as_posix(), level=logging.INFO, stdout=args.stdout)

    if args.platform == "gpu":
        if platform.system() == "Darwin":
            device = torch.device("mps")
        else:
            device = torch.device("cuda")
    else:
        device = torch.device("cpu")
        logging.info("Running on CPU / number of threads is set to 1")
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        torch.set_num_threads(1)
        logging.info("Number of threads: {}".format(torch.get_num_threads()))
        logging.info("Please specify the number of threads with OMP_NUM_THREADS and MKL_NUM_THREADS")

    logging.info(f"device = {device}")
    print_args_info(args)

    optim_name = args.optimizer
    iter = args.num_iter
    dtype = get_dtype(args.dtype)

    if args.loss == "qsmel":
        if "1D" in args.model:
            L = [args.length1]
            loss_dir = f"{params['lt']}_{L[0]}{args.loss}"
        elif "2D" in args.model:
            L = [args.length1, args.length2]
            loss_dir = f"{params['lt']}_{L[0]}x{L[1]}{args.loss}"
        else:
            raise ValueError("Invalid model name. Must contain 1D or 2D")
        H, _, _ = get_model(args.model, params, L)
        h_list = [H]
        loss = rms_torch.SystemQuasiEnergyLoss(h_list, device=device)
    elif args.loss == "stoq":
        h_list, _, _ = get_model(args.model, params)
        loss = rms_torch.SystemStoquastic(h_list, device=device)
        loss_dir = f"{params['lt']}_{args.loss}"
    elif args.loss == "mel":
        h_list, _, _ = get_model(args.model, params)
        loss = rms_torch.MinimumEnergyLoss(torch.tensor(h_list, dtype = dtype, device = device), device=device, decay=epochs/10, dtype=dtype)
        loss_dir = f"{params['lt']}_{args.loss}"
    elif args.loss == "none":
        h_list, _, _ = get_model(args.model, params)
        loss = rms_torch.MinimumEnergyLoss(torch.tensor(h_list), device=device)
        loss_dir = f"{params['lt']}_{args.loss}"
        optim_name = "none"
        logging.info(
            "Loss : None is specified. " +
            "Loss function is automatically set to mel and no optimization will be performed. " +
            "Iteration set to 0")
        iter = 0
    else:
        raise ValueError("Invalid loss function")

    seed_list = [np.random.randint(0, 1000000) for i in range(iter)]
    # seed_list = range(iter, 2*iter)
    # use system hamiltonian for qsmel

    local_h_list, sps, model_name = get_model(args.model, params)
    custom_dir = Path("out/tensorboard")
    custom_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists
    setting_name = f"lr_{args.learning_rate}_epoch_{epochs}"
    h_path = Path("array/torch") / model_name
    h_path.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists
    h_path = h_path / loss_dir
    u_path = h_path / optim_name / setting_name
    u_path.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists

    logging.info(f"Unitary will be saved to {u_path.resolve()}")
    logging.info(f"Hamiltonian saved to {h_path}/H/")

    # save local hamiltonian
    save_npy(h_path / "H", [-np.array(h) for h in local_h_list])  # minus for - beta * H
    # if args.loss == "none":
    #     identity = torch.eye(h_list[0].shape[1], dtype=torch.float64)
    #     initial_loss = loss(identity).item()
    #     logging.info(f"initial loss = {initial_loss}")
    #     exit(0)

    # save global hamiltonian

    optimizer_func: type[torch.optim.Optimizer] = rms_torch.LION
    if args.optimizer == "LION":
        optimizer_func = rms_torch.LION
        learning_params = dict(
            lr=args.learning_rate,
        )
    elif args.optimizer == "Adam":
        optimizer_func = rms_torch.Adam
        learning_params = dict(
            lr=args.learning_rate,
            betas=(0.3, 0.5),
        )
    else:
        raise ValueError("Invalid optimizer")

    model = rms_torch.UnitaryRiemann(
        h_list[0].shape[1], sps, device=device).to(device)
    model.reset_params(torch.eye(sps))

    initial_loss = loss(model()).item()
    logging.info(f"initial loss = {initial_loss}")

    best_loss = initial_loss
    best_us = [
        np.eye(sps, dtype=np.float64),
    ]
    num_print = 10
    for i, seed in enumerate(seed_list):
        logging.info(f"------ Start iteration {i + 1} ------")
        start = time.time()
        torch.manual_seed(seed)
        np.random.seed(seed)
        loss.initializer()
        local_best_loss = 1e10
        local_best_us = []
        model.reset_params()
        optimizer = optimizer_func(model.parameters(), **learning_params)
        logging.info(f"seed = {seed}")

        for t in range(epochs):
            optimizer.zero_grad()
            output = model()
            loss_val = loss(output)
            loss_val_item = loss_val.item()
            if loss_val_item < local_best_loss:
                with torch.no_grad():
                    local_best_loss = np.copy(loss_val_item)
                    local_best_us = [np.copy(p.data.detach().cpu().numpy())
                                     for p in model.parameters()]

            loss_val.backward()
            model.update_riemannian_gradient()
            if (t+1) % (epochs // num_print) == 0 or t == 0:
                logging.info(
                    f"I: {i + 1}/{len(seed_list)} : Epoch: {t+1}/{epochs}, Loss: {loss_val.item()}", )
            optimizer.step()

        if local_best_loss < best_loss:
            best_loss = np.copy(local_best_loss)
            best_us = [np.copy(u) for u in local_best_us]

        time_elapsed = time.time() - start
        logging.info(
            f"""best loss at epoch {epochs}: {local_best_loss}, best loss so far: {best_loss},
            time elapsed: {time_elapsed:.4f} seconds"""
        )

        u_path_epoch = u_path / f"loss_{local_best_loss:.7f}/u"
        u_path_epoch.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists
        save_npy(u_path_epoch, local_best_us)

    model.reset_params(torch.from_numpy(best_us[0]))
    out = model()
    loss.initializer()
    if (abs(best_loss - loss(out).item())) > 1e-8:
        logging.error(
            """
            The best loss and the actual loss do not match.
            Something is wrong with the optimization.
            best_loss: {} and actual loss: {}
            """.format(best_loss, loss(out).item()))

    u_path_epoch = u_path / f"loss_{best_loss:.7f}/u"
    save_npy(u_path_epoch, best_us)

    logging.info(f"best loss: {best_loss} / initial loss: {initial_loss}")
    logging.info(f"best loss was saved to {u_path}/loss_{best_loss:.7f}/u")
    logging.info(f"hamiltonian was saved to {h_path}/H")

    if args.symoblic_link is not None:
        symb_path = args.symoblic_link
        logging.info(f"Link is {symb_path.resolve()}")
        if symb_path.exists():
            logging.info(f"Remove existing link {symb_path}")
            symb_path.unlink()
        logging.info(f"Create symbolic link to {h_path.resolve()}")
        os.symlink(h_path.resolve(), symb_path)

    # Save the some information to the hamiltonian directory
    # First read the existing info.txt file. If it exists, take best_loss and
    # initial_loss from there
    new_info = True
    info_path = h_path / "info.txt"
    logging.info(f"Save information to {info_path.resolve()}")
    if (info_path).exists():
        with open(info_path) as f:
            info = f.read()
            match = re.search(r"best loss: (.+?) \/ initial loss: (.+)", info)
            if match:
                prev_best_loss = float(match.group(1))
                prev_initial_loss = float(match.group(2))
                logging.info(
                    "Previous best loss: {} and previous initial loss: {}".format(
                        prev_best_loss, prev_initial_loss))
                if prev_initial_loss != initial_loss:
                    logging.warning(
                        f"initial loss in info.txt ({prev_initial_loss}) does not match with the current initial loss ({initial_loss})")
                    logging.warning("Delete old info.txt and create a new one")
                if prev_best_loss > best_loss:
                    logging.info(
                        f"New best loss is found. Update info.txt file. Loss was {prev_best_loss} and now {best_loss}")
                else:
                    new_info = False
                    logging.info(
                        f"Best loss did not change. Do not update info.txt file. Loss was {prev_best_loss} and now {best_loss}")
            else:
                logging.warning(
                    "Could not find best_loss and initial_loss from info.txt. Delete old info.txt and create a new one")
    else:
        logging.info("info.txt does not exist. Create a new one")

    if new_info:
        with open(info_path, "w") as f:
            best_unitary_path = u_path / f"loss_{best_loss:.7f}" / "u"
            hamiltonian_path = h_path / "H"
            f.write(f"best loss: {best_loss} / initial loss: {initial_loss}\n")
            f.write(f"best loss was saved to {best_unitary_path.resolve()}\n")
            f.write(f"hamiltonian was saved to {hamiltonian_path.resolve()}\n")
