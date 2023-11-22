import argparse
import random

models = [
    "KH",
    "HXYZ",
    "HXYZ2D",
    "FF1D",
    "FF2D",
]
# minimum energy solver, quasi energy solver
loss_val = ["mel", "none"]

parser = argparse.ArgumentParser(
    description="exact diagonalization of shastry_surtherland")
parser.add_argument("-m", "--model", help="model (model) Name",
                    required=True, choices=models)
parser.add_argument("-Jz", "--coupling_z", help="coupling constant (Jz)",
                    type=float, default=1)  # SxSx + SySy +
parser.add_argument("-Jx", "--coupling_x",
                    help="coupling constant (Jx)", type=float)
parser.add_argument("-Jy", "--coupling_y",
                    help="coupling constant (Jy)", type=float)
parser.add_argument("-hx", "--mag_x", help="magnetic field",
                    type=float, default=0)
parser.add_argument("-hz", "--mag_z", help="magnetic field",
                    type=float, default=0)
parser.add_argument("-T", "--temperature", help="temperature", type=float)
parser.add_argument("-M", "--num_iter",
                    help="# of iterations", type=int, default=10)
parser.add_argument("-r", "--seed", help="random seed",
                    type=int, default=random.randint(0, 1000000))
parser.add_argument("--sps", help="sps", type=int, default=3)
parser.add_argument("-lr", "--learning_rate",
                    help="learning rate", type=float, default=0.01)
parser.add_argument(
    "--print", help="print out the result", action="store_true")
parser.add_argument("-e", "--epochs", help="epoch", type=int, default=100)
parser.add_argument(
    "-lt",
    "--lattice_type",
    help="algorithm determine local unitary matrix",
    default="original",
)
parser.add_argument(
    "-loss",
    "--loss",
    help="loss_methods",
    choices=loss_val,
    default="mel",
    nargs="?",
    const="all",
)

parser.add_argument(
    "-o",
    "--optimizer",
    help="optimizer",
    choices=["LION", "Adam"],
    default="LION",
    nargs="?",
    const="all",
)

parser.add_argument(
    "-p",
    "--platform",
    help="cput / gpu",
    choices=["cpu", "gpu"],
    nargs="?",
    const="all",
    default="cpu",
)

parser.add_argument(
    "-n",
    "--num_threads",
    help="number of threads",
    type=int,
    default=1,
)


args = parser.parse_args()

params = dict(
    Jx=args.coupling_x if args.coupling_x is not None else args.coupling_z,
    Jy=args.coupling_y if args.coupling_y is not None else args.coupling_z,
    Jz=args.coupling_z,
    hx=args.mag_x,
    hz=args.mag_z,
    lt=args.lattice_type,
    sps=args.sps,
    seed=args.seed, # random seed to generate hamiltonian
)


args_str = "args: {}".format(args)
hash_str = str(hash(args_str))
