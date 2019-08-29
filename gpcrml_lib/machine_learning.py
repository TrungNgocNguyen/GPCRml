from data import GetDescriptors
import pandas as pd
import pickle
import numpy as np
import seaborn as sns


class DecisionTreeModel(GetDescriptors):

    dtc_model = 'dtc_3_features.sav'

    def __init__(self, getdescriptors):
        self.dihedrals_phi = getdescriptors.dihedrals_phi
        self.dihedrals_psi = getdescriptors.dihedrals_psi

        self.phipsi_concat = None
        self.prediction = None

    def predict(self):
        self.phipsi_concat = pd.concat([self.dihedrals_psi, self.dihedrals_phi], axis=1)

        dtc = pickle.load(open(self.dtc_model, 'rb'))
        self.prediction = dtc.predict(self.phipsi_concat)

        unique, counts = np.unique(self.prediction, return_counts=True)
        return dict(zip(unique, counts))

    def plot_results(self):
        sns.set(style="darkgrid")
        values = self.prediction
        sns.countplot(x=values, linewidth=0.5)