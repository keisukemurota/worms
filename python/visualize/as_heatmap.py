"""plot average sign and loss as heatmap."""
import numpy as np  # noqa
# import os  # noqa
# import matplotlib.gridspec as gridspec  # noqa
import matplotlib.pyplot as plt  # noqa
import pandas as pd  # noqa
import itertools  # noqa
# import pandas  # noqa

import sys  # noqa
from pathlib import Path  # noqa
PYTHON_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, PYTHON_DIR.resolve().as_posix())
print(sys.path)
from rmsKit import utils  # noqa
logger = utils.get_logger("log.log", level="INFO", stdout=True)
parser = utils.parser.get_parser()
args, _, _ = utils.parser.get_params_parser(parser)

IMAGE_PATH = PYTHON_DIR / "visualize" / "image"
WORM_RESULT_PATH = PYTHON_DIR / "rmsKit" / "array" / "quetta"

model_name = args.model
image_model_dir = IMAGE_PATH / model_name
worm_result_path = WORM_RESULT_PATH / (model_name+"_loc")

N = 10**6
BETA_THRES = 20

if not IMAGE_PATH.exists():
    raise FileNotFoundError("{} does not exist.".format(IMAGE_PATH.resolve()))
if not WORM_RESULT_PATH.exists():
    raise FileNotFoundError("{} does not exist.".format(WORM_RESULT_PATH.resolve()))
if not worm_result_path.exists():
    raise FileNotFoundError("{} does not exist.".format(worm_result_path.resolve()))

summary_files = utils.path.find_summary_files(worm_result_path)
df = utils.path.get_df_from_summary_files(summary_files, N)
df = df[df.sweeps == N]
df = df[df["T"] >= 1 / BETA_THRES]
logger.info("temeprature simulated: {}".format(np.sort(df["T"].unique())))
logger.info("L simulated: {}".format(np.sort(df.n_sites.unique())))

params_df = utils.param_dict_normalize(df['ham_path'].apply(utils.extract_parameters_from_path))
df = pd.concat([df, params_df], axis=1)
df = df.rename(columns={"T": "temperature"})
print(df.columns)
print(len(df))

image_model_dir.mkdir(parents=False, exist_ok=True)

# Define the heatmap plotting function


# Define the heatmap plotting function
def plot_heatmap(df, fixed_params, heatmap_params, model_name, image_model_dir):
    """Plot the heatmap of the average sign and loss as a function of two parameters."""
    print(fixed_params)
    print(df.temperature.unique(), df.n_sites.unique())
    print(df.columns)
    df_t = df[(df.temperature == 2) & (df.n_sites == 25) & (df["Jy"] == -1.0) & (df["hx"] == 0.0)]
    df_t.to_csv(image_model_dir / "filtered_data.csv", index=False)
    # exit()
    for fixed_values in itertools.product(*fixed_params.values()):
        filtered_df = df.copy()
        figure_name_parts = []
        # print(fixed_values)
        for key, value in zip(fixed_params.keys(), fixed_values):
            filtered_df = filtered_df[filtered_df[key] == value]
            # print(key, value, len(filtered_df))
            figure_name_parts.append(f"{key}_{value}")
            print(f"{key} : {value}")
            print("after filtering : len for fixed values : ", len(filtered_df))
        


        x_param, y_param = heatmap_params
        x_values = np.sort(filtered_df[x_param].unique())
        y_values = np.sort(filtered_df[y_param].unique())
        x, y = np.meshgrid(x_values, y_values)

        zs = {
            "AS (optimized)": [],
            "AS (original)": [],
            "$log(\\eta) / \\beta$ (optimized)": [],
            "$log(\\eta) / \\beta$ (original)": [],
            "Loss (optimized)": [], "Loss (original)": [],
        }

        # Process data for heatmap
        # print(x, y)
        for Jx, Jy in zip(x.reshape(-1), y.reshape(-1)):
            df_plot = filtered_df[(filtered_df[x_param] == Jx) & (filtered_df[y_param] == Jy)]
            # df_u = df_plot[~df_plot.loss.isna()]
            df_u = df_plot

            if len(df_u) == 0:
                logger.info(
                    f"Optimized system data was not found for {x_param}={Jx}, {y_param}={Jy}")
                au = np.nan
                au_err = np.nan
                loss = np.nan
                init_loss = np.nan
            else:
                idx = np.argmin(df_u.loss.values)
                loss = df_u.loss.values[idx]
                init_loss = df_u.init_loss.values[idx]
                au = df_u["as"].values[idx]
                au_err = df_u["as_error"].values[idx]

                if np.abs(au) < np.abs(au_err) * 5:
                    au = np.nan
                    au_err = np.nan
                
                au_err *= np.sqrt(N)



            df_h = df_plot[df_plot.loss.isna()]
            if len(df_h) == 0:
                # logger.info(
                #     f"Original System data was not found for {x_param}={Jx}, {y_param}={Jy}")
                ah = np.nan
                ah_err = np.nan
            else:
                ah = df_h["as"].min()
                ah_err = df_h["as_error"].min() * np.sqrt(N)

            au = np.maximum(au, 1e-5)
            ah = np.maximum(ah, 1e-5)
            nu = -np.log(au) * filtered_df.temperature.values[0]
            nh = -np.log(ah) * filtered_df.temperature.values[0]

            # print(f"temperature : {filtered_df.temperature.values[0]}")
            zs["AS (optimized)"].append(au)
            zs["AS (original)"].append(ah)
            zs["$log(\\eta) / \\beta$ (optimized)"].append(nu)
            zs["$log(\\eta) / \\beta$ (original)"].append(nh)
            zs["Loss (optimized)"].append(loss)
            zs["Loss (original)"].append(init_loss)

        zs["AS (optimized)"] = np.array(zs["AS (optimized)"])
        zs["AS (original)"] = np.array(zs["AS (original)"])
        zs["Loss (optimized)"] = np.array(zs["Loss (optimized)"])
        zs["Loss (original)"] = np.array(zs["Loss (original)"])
        zs["$log(\\eta) / \\beta$ (optimized)"] = np.array(zs["$log(\\eta) / \\beta$ (optimized)"])
        zs["$log(\\eta) / \\beta$ (original)"] = np.array(zs["$log(\\eta) / \\beta$ (original)"])
        # print("loss optimized : ",zs["Loss (optimized)"])
        # Replace nan with largest value
        zs["AS (optimized)"][np.isnan(zs["AS (optimized)"])] = np.nanmin(
            zs["AS (optimized)"])
        zs["AS (original)"][np.isnan(zs["AS (original)"])] = np.nanmin(
            zs["AS (original)"])
        zs["Loss (optimized)"][np.isnan(zs["Loss (optimized)"])] = np.nanmin(
            zs["Loss (optimized)"])
        zs["Loss (original)"][np.isnan(zs["Loss (original)"])] = np.nanmin(
            zs["Loss (original)"])
        zs["$log(\\eta) / \\beta$ (optimized)"][np.isnan(zs["$log(\\eta) / \\beta$ (optimized)"])] = np.nanmax(
            zs["$log(\\eta) / \\beta$ (optimized)"])
        zs["$log(\\eta) / \\beta$ (original)"][np.isnan(zs["$log(\\eta) / \\beta$ (original)"])] = np.nanmax(
            zs["$log(\\eta) / \\beta$ (original)"])
        
        zs["$log(\\eta) / \\beta$ (optimized)"][np.isnan(zs["$log(\\eta) / \\beta$ (optimized)"])] = 100
        zs["$log(\\eta) / \\beta$ (original)"][np.isnan(zs["$log(\\eta) / \\beta$ (original)"])] = 100

        # print(zs["$log(\\eta) / \\beta$ (optimized)"])
        # print(zs["$log(\\eta) / \\beta$ (original)"])

        zs.pop("AS (optimized)")
        zs.pop("AS (original)")
        # print(zs["Loss (optimized)"])

        # Plotting
        fig, ax = plt.subplots(2, 2, figsize=(13, 13))
        try:
            max_loss = max(np.max(zs["Loss (optimized)"]), np.max(zs["Loss (original)"]))
        except BaseException as e:
            raise e

        max_eta = max(np.max(zs["$log(\\eta) / \\beta$ (optimized)"]), np.max(zs["$log(\\eta) / \\beta$ (original)"]))
        max_eta = min(max_eta, 3)
        max_as = 1
        if np.abs(max_loss) < 1e-5:
            max_loss = 1

        for i, (key, z) in enumerate(zs.items()):
            Z = np.array(z).reshape(x.shape)
            if "eta" in key:
                vmin, vmax = (0, 1.5)
            elif "Loss" in key:
                vmin, vmax = (0, max_loss)
            else:
                vmin, vmax = (0, max_as)

            # Choose colormap based on whether "Loss" is in the key
            if "Loss" in key:
                colormap = 'Reds'  # Colormap for loss plots
            elif "eta" in key:
                colormap = 'RdPu'  # Colormap for eta plots
            else:
                colormap = 'RdPu_r'  # Colormap for other plots
            if "eeeee" in key:
                c = ax[i % 2, i // 2].imshow(Z, cmap=colormap, aspect="auto",
                                            extent=[x.min(), x.max(), y.min(), y.max()],
                                            origin='lower', interpolation='none')
            else:
                c = ax[i % 2, i // 2].imshow(Z, cmap=colormap, vmin=vmin, vmax=vmax, aspect="auto",
                                            extent=[x.min(), x.max(), y.min(), y.max()],
                                            origin='lower', interpolation='none')
            ax[i % 2, i // 2].set_title(key, fontsize=25)
            ax[i % 2, i // 2].set_xlabel(x_param, fontsize=20)
            ax[i % 2, i // 2].set_ylabel(y_param, fontsize=20)
            # Adjust font size for tick labels
            ax[i % 2, i // 2].tick_params(axis='both', which='major', labelsize=15)
            cbar = fig.colorbar(c, ax=ax[i % 2, i // 2], fraction=0.06, pad=0.04)
            if "Loss" in key or "eta" in key:
                pass
            else:
                cbar.ax.invert_yaxis()

        plt.tight_layout()
        figure_name = f"AsHeatmap_{'_'.join(figure_name_parts)}.png"
        # figure_name = f"AsHeatmap_{'_'.join(figure_name_parts)}.pdf"
        figure_path = image_model_dir / figure_name
        # plt.savefig(figure_path, bbox_inches='tight', format="pdf")
        plt.savefig(figure_path, bbox_inches='tight')
        logger.info(f"Figure saved to {figure_path}")


if model_name == "HXYZ1D":
    fixed_params_HXYZ1D = {
        "temperature": np.sort(
            df.temperature.unique()), "n_sites": np.sort(
            df.n_sites.unique()), "hx": np.sort(
                df.hx.unique())}
    plot_heatmap(df, fixed_params_HXYZ1D, ('Jx', 'Jy'), model_name, image_model_dir)


if model_name == "HXYZ2D":
    fixed_params_HXYZ2D = {
        "temperature": np.sort(
            df.temperature.unique()), "n_sites": np.sort(
            df.n_sites.unique()), "hx": np.sort(
                df.hx.unique())}
    plot_heatmap(df, fixed_params_HXYZ2D, ('Jx', 'Jy'), model_name, image_model_dir)

elif model_name == "MG1D":
    fixed_params_MG1D = {
        "temperature": np.sort(df.temperature.unique()),
        "n_sites": np.sort(df.n_sites.unique()),
    }
    plot_heatmap(df, fixed_params_MG1D, ('J2', 'J3'), model_name, image_model_dir)

elif model_name == "BLBQ1D":
    fixed_params_BLBQ1D = {
        "temperature": np.sort(df.temperature.unique()),
        "n_sites": np.sort(df.n_sites.unique()),
        "J0": [1],
        "bc": ["obc"],
    }
    plot_heatmap(df, fixed_params_BLBQ1D, ('J1', 'hx'), model_name, image_model_dir)

elif model_name == "SS2D":
    fixed_params_MG1D = {
        "temperature": [0.1, 0.25, 1],
        "n_sites": [16,36],
        "J0": [1],
        "loss_func": ["-1_none", "1_mel"]
    }
    plot_heatmap(df, fixed_params_MG1D, ('J1', 'J2'), model_name, image_model_dir)

elif model_name == "KH2D":
    fixed_params_KH2D = {
        "temperature": [2, 4],
        "n_sites": [16, 25],
        "loss_func": ["3_mel"],
        "hx" : [0.0, 0.5],
    }
    plot_heatmap(df, fixed_params_KH2D, ('Jx', 'Jy'), model_name, image_model_dir)
