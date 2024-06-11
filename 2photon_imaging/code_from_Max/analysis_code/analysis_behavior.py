from raw_preprocess import data_fetcher, peak_extractor
import os
import os.path as path
from visualize_trials_comparison import trials_plot
from visualize_video import tmp_session_videos_interface
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde, ttest_ind
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from raw_preprocess import behavior_extractor, auc_extractor
from utils import simple_beeswarm2

color_dict = {
    "PPH_M_1R": "red",
    "PPH_M_NC": "blue",
}


def movement_vs_whisking(all_sessions):
    fig, ax = plt.subplots(1, 2, sharex='all', sharey='all')
    mice_1r_x, mice_1r_y = [], []
    mice_nc_x, mice_nc_y = [], []
    for sessions in all_sessions:
        for session in sessions:
            movement_trace = gaussian_filter(-session.movement * 10, sigma=1)
            whisking_trace = gaussian_filter(session.whisking * 100, sigma=1)

            if session.exp_mice == "PPH_M_1R":
                mice_1r_x.append(movement_trace)
                mice_1r_y.append(whisking_trace)
            else:
                mice_nc_x.append(movement_trace)
                mice_nc_y.append(whisking_trace)
    mice_1r_x, mice_1r_y = np.concatenate(mice_1r_x), np.concatenate(mice_1r_y)
    mice_nc_x, mice_nc_y = np.concatenate(mice_nc_x), np.concatenate(mice_nc_y)
    xy = np.vstack([mice_1r_x, mice_1r_y])
    z = gaussian_kde(xy)(xy)
    ax[0].scatter(mice_1r_x, mice_1r_y, marker='.', c=z, s=5, alpha=0.6)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax[0])
    cbar.ax.set_ylabel('Density')
    # ax[0].hist2d(mice_1r_x, mice_1r_y, bins=(100, 100), cmap=plt.cm.jet)
    xy = np.vstack([mice_nc_x, mice_nc_y])
    z = gaussian_kde(xy)(xy)
    ax[1].scatter(mice_nc_x, mice_nc_y, marker='.', c=z, s=5, alpha=0.6)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax[1])
    cbar.ax.set_ylabel('Density')
    # ax[1].hist2d(mice_nc_x, mice_nc_y, bins=(100, 100), cmap=plt.cm.jet)
    for i in range(2):
        ax[i].spines[['right', 'top']].set_visible(False)
        ax[i].set_xlabel("Velocity [cm/s]")
        ax[i].set_ylabel("Whisking Intensity [%]")
    ax[0].set_title("PPH_M_1R")
    ax[1].set_title("PPH_M_NC")
    plt.show()
    exit()


def movement_vs_whisking_per_trial(all_trials):
    fig, ax = plt.subplots(1, 2, sharex='all', sharey='all')
    mice_1r_x, mice_1r_y = [], []
    mice_nc_x, mice_nc_y = [], []
    for trials in all_trials:
        for trial in trials:
            if not trial.drop:
                movement_auc = np.sum(np.abs(trial.movement[20:30]))
                whisking_auc = np.sum(np.abs(trial.whisking[20:30]))
                if trial.exp_mice == "PPH_M_1R":
                    mice_1r_x.append(movement_auc)
                    mice_1r_y.append(whisking_auc)
                else:
                    mice_nc_x.append(movement_auc)
                    mice_nc_y.append(whisking_auc)
    mice_1r_x, mice_1r_y = np.array(mice_1r_x), np.array(mice_1r_y)
    mice_nc_x, mice_nc_y = np.array(mice_nc_x), np.array(mice_nc_y)
    xy = np.vstack([mice_1r_x, mice_1r_y])
    z = gaussian_kde(xy)(xy)
    ax[0].scatter(mice_1r_x, mice_1r_y, marker='.', c=z, s=20, alpha=0.6)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax[0])
    cbar.ax.set_ylabel('Density')
    # ax[0].hist2d(mice_1r_x, mice_1r_y, bins=(100, 100), cmap=plt.cm.jet)
    xy = np.vstack([mice_nc_x, mice_nc_y])
    z = gaussian_kde(xy)(xy)
    ax[1].scatter(mice_nc_x, mice_nc_y, marker='.', c=z, s=20, alpha=0.6)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax[1])
    cbar.ax.set_ylabel('Density')
    # ax[1].hist2d(mice_nc_x, mice_nc_y, bins=(100, 100), cmap=plt.cm.jet)
    for i in range(2):
        ax[i].spines[['right', 'top']].set_visible(False)
        ax[i].set_xlabel("Velocity AUC [cm]")
        ax[i].set_ylabel("Whisking Intensity AUC [Frames]")
        ax[i].axhline(y=4, lw=1, ls='--', color='black')
        ax[i].axhline(y=9.5, xmin=0.1, lw=1, ls='--', color='black')
        ax[i].axvline(x=1, lw=1, ls='--', color='black')
    ax[0].set_title("PPH_M_1R")
    ax[1].set_title("PPH_M_NC")
    plt.show()
    exit()


def auc_ttest(all_stationary, all_whi_only, all_whi_while, all_intense):
    name_list = [
        "stationary",
        "whisking only",
        "whisking while running",
        "intense whisking",
    ]
    fig, ax = plt.subplots(1, 2, sharex='all', sharey='all')
    a_nc, a_1r, b_nc, b_1r = [], [], [], []
    for stationary, whisking_only, whisking_while_running, intense_whisking in \
            zip(all_stationary, all_whi_only, all_whi_while, all_intense):
        if len(intense_whisking) > 0:
            auc_list = [
                stationary.extract("df_f_auc"),
                # whisking_only.extract("df_f_auc"),
                np.concatenate([
                    whisking_while_running.extract("df_f_auc"),
                    intense_whisking.extract("df_f_auc")], axis=0)
            ]
        else:
            auc_list = [
                stationary.extract("df_f_auc"),
                # whisking_only.extract("df_f_auc"),
                whisking_while_running.extract("df_f_auc"),
            ]
        if stationary[0].exp_mice == "PPH_M_1R":
            a_1r.append(auc_list[0].reshape(-1))
            b_1r.append(auc_list[1].reshape(-1))
        else:
            a_nc.append(auc_list[0].reshape(-1))
            b_nc.append(auc_list[1].reshape(-1))
    a_nc, b_nc = np.concatenate(a_nc), np.concatenate(b_nc)
    a_1r, b_1r = np.concatenate(a_1r), np.concatenate(b_1r)
    p = ttest_ind(a_nc, b_nc).pvalue
    print(p)
    ax[1].scatter(simple_beeswarm2(a_nc, width=0.25) + 0, a_nc, s=5, facecolors='none', edgecolors='black', alpha=0.7)
    ax[1].errorbar(0 - 0.5, np.mean(a_nc), yerr=np.std(a_nc),
                   fmt='o', color='black', alpha=0.7, markersize=4, capsize=7)
    ax[1].scatter(simple_beeswarm2(b_nc, width=0.25) + 2, b_nc, s=5, facecolors='none', edgecolors='black', alpha=0.7)
    ax[1].errorbar(2 - 0.5, np.mean(b_nc), yerr=np.std(b_nc),
                   fmt='o', color='black', alpha=0.7, markersize=4, capsize=7)
    ax[1].set_ylabel(r"Trial AUC [df/f_0 * frames]")
    ax[1].set_xlim(-1, 2.5)

    p = ttest_ind(a_1r, b_1r).pvalue
    ax[0].scatter(simple_beeswarm2(a_1r, width=0.25) + 0, a_1r, s=5, facecolors='none', edgecolors='black', alpha=0.7)
    ax[0].errorbar(0 - 0.5, np.mean(a_1r), yerr=np.std(a_1r),
                   fmt='o', color='black', alpha=0.7, markersize=4, capsize=7)
    ax[0].scatter(simple_beeswarm2(b_1r, width=0.25) + 2, b_1r, s=5, facecolors='none', edgecolors='black', alpha=0.7)
    ax[0].errorbar(2 - 0.5, np.mean(b_1r), yerr=np.std(b_1r),
                   fmt='o', color='black', alpha=0.7, markersize=4, capsize=7)
    ax[0].set_ylabel(r"Trial AUC [df/f_0 * frames]")
    ax[0].set_xlim(-1, 2.5)
    for i in range(2):
        ax[i].spines[['right', 'top']].set_visible(False)
        ax[i].set_xticks([0, 2], ["stationary", f"whisking while running\n+\nintense whisking"])
    print(p)
    plt.show()


def main(task_id, analysis_list, data_path, save_dir):
    all_trials, all_sessions = [], []
    all_stationary, all_whi_only, all_whi_while, all_intense = [], [], [], []
    for mice_id, num_session in analysis_list:
        trials, sessions = data_fetcher(path.join(data_path, task_id, mice_id), num_session,
                                        path.join(save_dir, task_id, mice_id), elaborate=False, )
        all_trials.append(trials)
        all_sessions.append(sessions)
        tmp_session_videos_interface(sessions, path.join(data_path, task_id, mice_id),
                                     save_dir=path.join(save_dir, task_id, mice_id, "videos"))
        continue
        trials.append(auc_extractor)
        trials.append(behavior_extractor)
        stationary = trials.filter(drop=False, behavior_type='stationary')
        whisking_only = trials.filter(drop=False, behavior_type='whisking only')
        whisking_while_running = trials.filter(drop=False, behavior_type='whisking while running')
        intense_whisking = trials.filter(drop=False, behavior_type='intense whisking')

        all_stationary.append(stationary)
        all_whi_only.append(whisking_only)
        all_whi_while.append(whisking_while_running)
        all_intense.append(intense_whisking)

        # trials_plot([stationary, whisking_only, whisking_while_running, intense_whisking],
        #             [mice_id, ],
        #             ["stationary", "whisking only", "whisking while running", "intense whisking"],
        #             path.join(save_dir, task_id, mice_id, "behavior_compare.jpg"))
    # movement_vs_whisking(all_sessions)
    # movement_vs_whisking_per_trial(all_trials)
    # auc_ttest(all_stationary, all_whi_only, all_whi_while, all_intense)


if __name__ == "__main__":
    main("BehaviorAnalysis",
         (
             # (path.join("PPH_M_1R", "FOV1"), 2),
             # (path.join("PPH_M_1R", "FOV2"), 2),
             # (path.join("PPH_M_1R", "FOV3"), 2),
             (path.join("PPH_M_1R", "FOV4"), 2),
             # (path.join("PPH_M_NC", "FOV1"), 2),
             # (path.join("PPH_M_NC", "FOV2"), 2),
             # (path.join("PPH_M_NC", "FOV3"), 2),
             # (path.join("PPH_M_NC", "FOV4"), 2),
         ),
         path.join(path.dirname(__file__), "data"),
         path.join(path.dirname(__file__), "imgs", 'behaviors')
         )
