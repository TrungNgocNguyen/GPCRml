# python standard packages
import pickle

# external packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


class DecisionTreePredictions:
    def __init__(self, df_phi, df_psi):
        self.df_phi = df_phi
        self.df_psi = df_psi
        self.dtc = pickle.load(open("gpcrml_lib/dtc_3_features.sav", 'rb'))

        self.df_predictions = None

    def predict(self):
        phipsi_concat = pd.concat([self.df_psi, self.df_phi], axis=1)
        phipsi_concat = phipsi_concat[["3x52_phi", "7x54_phi", "6x47_psi"]]

        prediction = self.dtc.predict(phipsi_concat)

        phipsi_concat['Active(1)/Inactive(0)'] = prediction
        active_inactive_prediction = phipsi_concat.drop(["3x52_phi", "7x54_phi", "6x47_psi"], axis=1)
        self.df_predictions = active_inactive_prediction

        unique, counts = np.unique(prediction, return_counts=True)
        total_states = len(prediction)
        if len(counts) > 1:
            active_percent = counts[1]/total_states*100
            inactive_percent = counts[0]/total_states*100
            print("Predicted as Active State Conformations: {} ({}%)".format(counts[1], active_percent))
            print("Predicted as Inactive State Conformations: {} ({}%)".format(counts[0], inactive_percent))
            print("Total amount: {}".format(total_states))
        elif unique == 1:
            print("All conformations ({}) are predicted as Active State".format(counts[0]))
        else:
            print("All conformations ({}) are predicted as Inactive State".format(counts[0]))
        self.df_predictions = active_inactive_prediction
        return active_inactive_prediction

    def plot_predictions(self):
        sns.set(style="darkgrid")

        if self.df_predictions is None:
            print("No prediction data found. Use 'instance.predict()' to calculate predictions first.")
        else:
            values = self.df_predictions['Active(1)/Inactive(0)']
            sns.countplot(x=values, linewidth=0.5)
            plt.show()

    def plot_dihedrals(self):
        phipsi_concat = pd.concat([self.df_psi, self.df_phi], axis=1)
        bwn_columns = ["3x52_phi", "7x54_phi", "6x47_psi"]
        values = phipsi_concat[bwn_columns]

        sns.set(style="whitegrid")
        data = pd.DataFrame(values, columns=bwn_columns)
        data = data.rolling(1).mean()
        sns.lineplot(data=data, palette="tab10", linewidth=0.5)

        plt.show()