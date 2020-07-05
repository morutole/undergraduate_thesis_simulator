import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_raw_csv():
    raw_simulate_df = pd.read_csv("sample.csv")
    search_word = "sat_position"
    raw_position_df = raw_simulate_df.loc[:,raw_simulate_df.columns.str.contains(search_word)]
    search_word = "sat_velocity_i" 
    raw_velocity_df = raw_simulate_df.loc[:,raw_simulate_df.columns.str.contains(search_word)]

    return raw_position_df.values, raw_velocity_df.values

def read_general_csv(file_name):
    value_df = pd.read_csv(file_name, header = None)
    value_nd = value_df.values
    position_nd = value_nd[:, 0:3]
    velocity_nd = value_nd[:, 3:6]

    return position_nd, velocity_nd

def extract_observed_correspondent_value(GPS_period, true_position_nd, true_velocity_nd):
    true_time_scale = np.linspace(0, simulation_time, true_position_nd.shape[0]).astype(int)
    ind = np.where(true_time_scale%GPS_period == 0)[0]
    true_observed_position_nd = true_position_nd[ind]
    true_observed_velocity_nd = true_velocity_nd[ind]

    return true_observed_position_nd, true_observed_velocity_nd

def get_label(index, kind):
    if index == 0:
        res = "x"
    elif index == 1:
        res = "y"
    elif index == 2:
        res = "z"
    
    if kind == "velocity":
        return "V" + res + "[m/s]"
    else:
        return res + "[m]"

raw_position_nd, raw_velocity_nd = read_raw_csv()
true_position_nd, true_velocity_nd = read_general_csv("true.csv")
observed_position_nd, observed_velocity_nd = read_general_csv("observed.csv")
estimate_position_nd, estimate_velocity_nd = read_general_csv("estimate.csv")

simulation_time = true_position_nd.shape[0] - 1 #開始時点まで含んでいるため
raw_time_scale = np.linspace(0, simulation_time, raw_position_nd.shape[0])
assert simulation_time%(observed_position_nd.shape[0] - 1) == 0, "GPS_periodがおかしい"
GPS_period = simulation_time/(observed_position_nd.shape[0] - 1)

true_observed_position_nd, true_observed_velocity_nd = extract_observed_correspondent_value(GPS_period, true_position_nd, true_velocity_nd)

observed_time_scale = np.linspace(0, simulation_time, observed_position_nd.shape[0])

def plot_absolute_graph():
    for i in range(3):
        plt.plot(raw_time_scale, raw_position_nd[:, i], label = "true")
        plt.plot(raw_time_scale, estimate_position_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "position"))
        plt.legend()
        plt.show()

    for i in range(3):
        plt.plot(raw_time_scale, raw_velocity_nd[:, i], label = "true")
        plt.plot(raw_time_scale, estimate_velocity_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "velocity"))
        plt.legend()
        plt.show()

def plot_difference_graph():
    for i in range(3):
        plt.scatter(observed_time_scale, observed_position_nd[:, i] - true_observed_position_nd[:, i], label = "observed")
        plt.scatter(raw_time_scale, estimate_position_nd[:, i] - raw_position_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "position"))
        plt.legend()
        plt.ylim(-2.0, 2.0)
        plt.show()

    for i in range(3):
        plt.scatter(observed_time_scale, observed_velocity_nd[:, i] - true_observed_velocity_nd[:, i], label = "observed")
        plt.scatter(raw_time_scale, estimate_velocity_nd[:, i] - raw_velocity_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "velocity"))
        plt.legend()
        plt.ylim(-0.2, 0.2)
        plt.show()

def plot_rms_graph():
    plt.scatter(observed_time_scale, np.linalg.norm(observed_position_nd - true_observed_position_nd, axis = 1), label = "observed")
    plt.scatter(raw_time_scale, np.linalg.norm(estimate_position_nd - raw_position_nd, axis = 1), label = "estimate")
    plt.xlabel("time[s]")
    plt.ylabel("position[m]")
    plt.legend()
    plt.show()

    plt.scatter(observed_time_scale, np.linalg.norm(observed_velocity_nd - true_observed_velocity_nd, axis = 1), label = "observed")
    plt.scatter(raw_time_scale, np.linalg.norm(estimate_velocity_nd - raw_velocity_nd, axis = 1), label = "estimate")
    plt.xlabel("time[s]")
    plt.ylabel("velocity[m/s]")
    plt.legend()
    plt.show()

plot_difference_graph()
plot_rms_graph()