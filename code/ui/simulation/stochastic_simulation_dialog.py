from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QFormLayout, QLineEdit, QPushButton

from models.simulation_settings import SimulationSettings
from simulation.gillespie_simulator import GillespieSimulator
from ui import common_widgets
from ui.gene_controller import GeneController


class StochasticSimulationDialog(QDialog):
    def _ok_button_clicked(self):
        time_text = self.time_field.text().strip()

        end_time = int(time_text) if time_text else 1

        species = list()
        for x in range(0, self.species_checkboxes.count()):
            checkbox = self.species_checkboxes.itemAt(x).widget()
            if checkbox.isChecked():
                species.append(checkbox.text())

        self.close()

        # Precision field does not apply to stochastic simulation!
        s = SimulationSettings(0, end_time, 0, [s.strip() for s in species])
        GillespieSimulator.visualise(GillespieSimulator.simulate(GeneController.get_instance().network, s), s)

    def __init__(self):
        super().__init__()
        main = QVBoxLayout()
        fields = QFormLayout()

        self.time_field = QLineEdit()
        fields.addRow(QLabel("Simulation time"), self.time_field)

        self.species_checkboxes = common_widgets.make_species_checkboxes_layout()

        self.ok_button = QPushButton("Ok")
        self.ok_button.clicked.connect(self._ok_button_clicked)

        main.addLayout(fields)
        main.addLayout(self.species_checkboxes)
        main.addWidget(self.ok_button)

        self.setLayout(main)
        self.setFixedHeight(140)
        self.setMinimumWidth(200)
        self.setWindowTitle("Simulation settings")
        self.exec_()
