from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QCheckBox, QGridLayout

from models.simulation_settings import SimulationSettings
from simulation.ode_simulator import OdeSimulator
from ui import common_widgets
from ui.gene_controller import GeneController


class DeterministicSimulationDialog(QDialog):

    def _ok_button_clicked(self):
        time_text = self.time_field.text().strip()
        sampling_text = self.sampling_field.text().strip()

        end_time = int(time_text) if time_text else 1

        species = list()
        for x in range(0, self.species_checkboxes.count()):
            checkbox = self.species_checkboxes.itemAt(x).widget()
            if checkbox.isChecked():
                species.append(checkbox.text())

        sampling_rate = int(sampling_text)

        self.close()
        s = SimulationSettings(0, end_time, sampling_rate, [s.strip() for s in species])
        o = OdeSimulator(GeneController.get_instance().network, s)
        o.visualise(o.simulate())

    def __init__(self):
        super().__init__()
        main = QVBoxLayout()
        fields = QFormLayout()

        self.time_field = QLineEdit()
        fields.addRow(QLabel("Simulation time"), self.time_field)

        self.sampling_field = QLineEdit()
        fields.addRow(QLabel("Sampling rate"), self.sampling_field)

        self.species_checkboxes = common_widgets.make_species_checkboxes_layout()

        self.ok_button = QPushButton("Ok")
        self.ok_button.clicked.connect(self._ok_button_clicked)

        main.addLayout(fields)
        main.addWidget(QLabel("Which species to show: "))
        main.addLayout(self.species_checkboxes)
        main.addWidget(self.ok_button)

        self.setMinimumWidth(200)
        self.setWindowTitle("Simulation settings")
        self.setLayout(main)
        self.exec_()
