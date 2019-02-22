from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout

from models.simulation_settings import SimulationSettings
from simulation.ode_simulator import OdeSimulator
from ui.gene_controller import GeneController


class DeterministicSimulationDialog(QDialog):

    def _ok_button_clicked(self):
        time_text = self.time_field.text().strip()
        species_text = self.species_field.text().strip()
        sampling_text = self.sampling_field.text().strip()

        end_time = int(time_text) if time_text else 1
        species = species_text.split(",")
        sampling_rate = int(sampling_text)

        s = SimulationSettings(0, end_time, sampling_rate, [s.strip() for s in species])
        o = OdeSimulator(GeneController.get_instance().network, s)
        o.visualise(o.simulate())

        self.close()

    def __init__(self):
        super().__init__()

        main = QVBoxLayout()

        fields = QFormLayout()

        self.time_field = QLineEdit()
        fields.addRow(QLabel("Simulation time"), self.time_field)

        self.sampling_field = QLineEdit()
        fields.addRow(QLabel("Sampling rate"), self.sampling_field)

        self.species_field = QLineEdit()
        fields.addRow(QLabel("Which species to show"), self.species_field)

        self.ok_button = QPushButton("Ok")
        self.ok_button.clicked.connect(self._ok_button_clicked)

        main.addWidget(QLabel("Species are comma separated"))
        main.addLayout(fields)
        main.addWidget(self.ok_button)

        self.setFixedHeight(180)
        self.setMinimumWidth(200)
        self.setWindowTitle("Simulation settings")
        self.setLayout(main)
        self.exec_()
