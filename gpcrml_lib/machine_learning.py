# python standard packages
import pickle

# external packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


class Predictions:
    r"""The Predictions class is a parent class which predicts the activation state (active/inactive) based on the used machine learning model.
    The Prediction class always requires a phi and a psi angle dataframe mapped to generic numbers of a GPCR structure.
    The extraction and mapping of phi and psi angles is provided via the Descriptors class in GPCRml/gpcrml_lib/data.py.
    Basic plotting functions are also provided in this class.

    :param df_phi:
        Pandas DataFrame of phi angles mapped to generic numbers

    :param df_psi:
        Pandas DataFrame of psi angles mapped to generic numbers

    :param model:
        machine learning model used for predictions

    :param df_predictions:
        prediction results stored in a Pandas DataFrame

    :param features:
        required features for predictions with machine learning model
    """

    def __init__(self, df_phi, df_psi):
        self.df_phi = df_phi
        self.df_psi = df_psi
        self.model = None
        self.df_predictions = None
        self.features = None

    def predict(self):
        r"""Predicts the GPCR activation state (active/inactive) using given model with required phi and psi angles.
        Results are stored in self.df_predictions.

        :returns:
            Pandas DataFrame with active/inactive prediction results
        """

        phipsi_concat = pd.concat([self.df_psi, self.df_phi], axis=1)
        phipsi_concat = phipsi_concat[self.features]

        prediction = self.model.predict(phipsi_concat)

        phipsi_concat['Active(1)/Inactive(0)'] = prediction

        active_inactive_prediction = phipsi_concat.drop(self.features, axis=1)
        self.df_predictions = active_inactive_prediction

        active_percent = None
        unique, counts = np.unique(prediction, return_counts=True)
        total_states = len(prediction)
        if len(counts) > 1:
            active_percent = counts[1]/total_states*100
            inactive_percent = counts[0]/total_states*100
            print("Predicted as Active State Conformations: {} ({}%)".format(counts[1], active_percent))
            print("Predicted as Inactive State Conformations: {} ({}%)".format(counts[0], inactive_percent))
            print("Total amount: {}".format(total_states))
        elif unique == 1:
            active_percent = 100
            print("All conformations ({}) are predicted as Active State".format(counts[0]))
        else:
            active_percent = 0
            print("All conformations ({}) are predicted as Inactive State".format(counts[0]))
        self.df_predictions = active_inactive_prediction
        return active_inactive_prediction

    def plot_predictions(self):
        r"""Creates a simple barplot of active/inactive predictions."""

        sns.set(style="darkgrid")

        if self.df_predictions is None:
            print("No prediction data found. Use 'instance.predict()' to calculate predictions first.")
        else:
            values = self.df_predictions['Active(1)/Inactive(0)']
            sns.countplot(x=values, linewidth=0.5)
            plt.show()

    def plot_dihedrals(self):
        r"""Creates a simple lineplot of used phi and/or psi angles used for predictions."""

        phipsi_concat = pd.concat([self.df_psi, self.df_phi], axis=1)
        bwn_columns = self.features
        values = phipsi_concat[bwn_columns]

        sns.set(style="whitegrid")
        data = pd.DataFrame(values, columns=bwn_columns)
        data = data.rolling(1).mean()
        sns.lineplot(data=data, palette="tab10", linewidth=0.5)

        plt.show()


class DTCPredictions(Predictions):
    r"""The DTCPredictions class is a child class inheriting all functionality from the Predictions parent class. It uses
    the DTC machine learning model to predict the activation states of GPCRs based on the phi and psi angles used by this model.
    """

    def __init__(self, df_phi, df_psi):
        super().__init__(df_phi, df_psi)
        self.model = pickle.load(open("/home/gocky/python/GPCRml/gpcrml_lib/DTC_ml_model.sav", 'rb'))
        self.features = ['1x56_psi', '2x39_psi', '4x51_psi', '4x55_psi', '6x44_psi', '7x47_psi', '7x48_psi', '7x52_phi', '7x53_phi', '7x54_phi']

class KNNPredictions(Predictions):
    r"""The KNNPredictions class is a child class inheriting all functionality from the Predictions parent class. It uses
    the KNN machine learning model to predict the activation states of GPCRs based on the phi and psi angles used by this model.
    """

    def __init__(self, df_phi, df_psi):
        super().__init__(df_phi, df_psi)
        self.model = pickle.load(open("/home/gocky/python/GPCRml/gpcrml_lib/KNN_ml_model.sav", 'rb'))
        self.features = ['1x56_psi', '2x39_psi', '4x51_psi', '4x55_psi', '6x44_psi', '7x47_psi', '7x48_psi', '7x52_phi', '7x53_phi', '7x54_phi']

class SVMPredictions(Predictions):
    r"""The SVMPredictions class is a child class inheriting all functionality from the Predictions parent class. It uses
    the SVM machine learning model to predict the activation states of GPCRs based on the phi and psi angles used by this model.
    """

    def __init__(self, df_phi, df_psi):
        super().__init__(df_phi, df_psi)
        self.model = pickle.load(open("/home/gocky/python/GPCRml/gpcrml_lib/SVM_ml_model.sav", 'rb'))
        self.features = ['1x56_psi', '2x39_psi', '4x51_psi', '4x55_psi', '6x44_psi', '7x47_psi', '7x48_psi', '7x52_phi', '7x53_phi', '7x54_phi']