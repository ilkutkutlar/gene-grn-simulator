from typing import List

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QComboBox, QTabWidget
from PyQt5.QtWidgets import QDialog, QLineEdit
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QListWidget, QPushButton
from PyQt5.QtWidgets import QWidget, QApplication, QMainWindow, QAction, QMessageBox, \
    QFileDialog

from input_output.sbml_parser import SbmlParser
from models.formulae import TranscriptionFormula, TranslationFormula, DegradationFormula, CustomFormula
from models.network import Network
from models.reaction import Reaction
from models.simulation_settings import SimulationSettings
from simulation.gillespie_simulator import GillespieSimulator
from simulation.ode_simulator import OdeSimulator
from ui.add_reaction_dialog2 import AddReactionDialog2
from ui.gene_controller import GeneController
from ui.reactions_tab import ReactionsTab
from ui.species_tab import SpeciesTab


def validate_species(species):
    if species:
        for s in species:
            if s not in GeneController.get_instance().network.species:
                return False
        return True
    else:
        # To ensure that empty sets are accepted
        return True


def show_error_message(message) -> bool:
    error_message = QMessageBox()
    error_message.setIcon(QMessageBox.Warning)
    error_message.setWindowTitle("Error")
    error_message.setStandardButtons(QMessageBox.Ok)
    error_message.setText(message)

    button = error_message.exec_()
    if button == QMessageBox.Ok:
        return True
    else:
        return False


class GeneWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout()

        tabs = QTabWidget()
        tabs.addTab(SpeciesTab(), "Species")
        tabs.addTab(ReactionsTab(), "Reactions")
        self._init_menubar()

        layout.addWidget(tabs)
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)
        self.setWindowTitle("Gene")
        self.show()

    def _init_menubar(self):
        self.menubar = self.menuBar()

        # File menu
        file = self.menubar.addMenu("File")
        open_file = QAction("Open SBML file", self)
        open_file.triggered.connect(self._handler_open_sbml_clicked)
        file.addAction(open_file)

        # Simulate menu
        simulate = self.menubar.addMenu("Simulate")

        deterministic = QAction("Deterministic (ODE)", self)
        stochastic = QAction("Stochastic (Gillespie Algorithm)", self)

        deterministic.triggered.connect(self._handler_deterministic_clicked)
        stochastic.triggered.connect(self._handler_stochastic_clicked)

        simulate.addAction(deterministic)
        simulate.addAction(stochastic)

    def _handler_open_sbml_clicked(self):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)
        # TODO: Why is the filter not working?
        # file_dialog.setFilter("SBML files (*.xml)")

        if file_dialog.exec_():
            # selectedFiles returns all files, we only want to open one
            filename: str = file_dialog.selectedFiles()[0]
            net: Network = SbmlParser.parse(filename)
            GeneController.get_instance().network = net

            # Display a nice message showing the contents of the file loaded
            message = QMessageBox()
            message.setText("SBML file has been opened and this network has been loaded: \n\n" + net.__str__())
            message.exec_()

            self._refresh_reactions_list(self.reactions_panel.m_list)
            self._refresh_species_list(self.species_panel.m_list)

    def _handler_deterministic_clicked(self):
        dia = QDialog()
        main = QVBoxLayout()

        info = QLabel("Species are comma separated")

        time_field = QLineEdit()
        time_field.setPlaceholderText("Simulation time")

        sampling_field = QLineEdit()
        sampling_field.setPlaceholderText("Sampling rate")

        species_field = QLineEdit()
        species_field.setPlaceholderText("Which species to show")

        ok_button = QPushButton("Ok")

        def dialog_ok_button_clicked_handler():
            time_text = time_field.text().strip()
            species_text = species_field.text().strip()
            sampling_text = sampling_field.text().strip()

            end_time = int(time_text) if time_text else 1
            species = species_text.split(",")
            sampling_rate = int(sampling_text)

            s = SimulationSettings(0, end_time, sampling_rate, [s.strip() for s in species])
            o = OdeSimulator(GeneController.get_instance().network, s)
            o.visualise(o.simulate())

            dia.close()

        ok_button.clicked.connect(dialog_ok_button_clicked_handler)

        main.addWidget(info)
        main.addWidget(time_field)
        main.addWidget(sampling_field)
        main.addWidget(species_field)
        main.addWidget(ok_button)
        dia.setFixedHeight(180)
        dia.setMinimumWidth(200)
        dia.setWindowTitle("Simulation settings")
        dia.setLayout(main)
        dia.exec_()

    def _handler_stochastic_clicked(self):
        dia = QDialog()
        main = QVBoxLayout()

        info = QLabel("Species are comma separated")

        time_field = QLineEdit()
        time_field.setPlaceholderText("Simulation time")

        species_field = QLineEdit()
        species_field.setPlaceholderText("Which species to show")

        ok_button = QPushButton("Ok")

        def dialog_ok_button_clicked_handler():
            time_text = time_field.text().strip()
            species_text = species_field.text().strip()

            end_time = int(time_text) if time_text else 1
            species = species_text.split(",")

            # Precision field does not apply to stochastic simulation!
            s = SimulationSettings(0, end_time, 0, [s.strip() for s in species])
            GillespieSimulator.visualise(GillespieSimulator.simulate(GeneController.get_instance().network, s), s)

            dia.close()

        ok_button.clicked.connect(dialog_ok_button_clicked_handler)

        main.addWidget(info)
        main.addWidget(time_field)
        main.addWidget(species_field)
        main.addWidget(ok_button)
        dia.setFixedHeight(140)
        dia.setMinimumWidth(200)
        dia.setWindowTitle("Simulation settings")
        dia.setLayout(main)
        dia.exec_()

    def _handler_remove_reactions_button(self):
        del GeneController.get_instance().network.reactions[self.reactions_panel.m_list.currentRow()]
        self._refresh_reactions_list(self.reactions_panel.m_list)

    def _handler_remove_species_button(self):
        species_id = self.species_panel.m_list.property("id" + str(self.species_panel.m_list.currentRow()))
        del GeneController.get_instance().network.species[species_id]
        self._refresh_species_list(self.species_panel.m_list)


app = QApplication([])
g = GeneWindow()
app.exec_()
