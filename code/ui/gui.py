from PyQt5.QtWidgets import QDialog, QLineEdit
from PyQt5.QtWidgets import QTabWidget
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QPushButton
from PyQt5.QtWidgets import QWidget, QApplication, QMainWindow, QAction, QMessageBox, \
    QFileDialog

from input_output.sbml_parser import SbmlParser
from models.network import Network
from models.simulation_settings import SimulationSettings
from simulation.gillespie_simulator import GillespieSimulator
from simulation.ode_simulator import OdeSimulator
from ui.deterministic_simulation_dialog import DeterministicSimulationDialog
from ui.gene_controller import GeneController
from ui.reactions_tab import ReactionsTab
from ui.reverse_engineering_modify_tab import ReverseEngineeringModifyTab
from ui.species_tab import SpeciesTab
from ui.stochastic_simulation_dialog import StochasticSimulationDialog


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

        self.species_tab = SpeciesTab()
        self.reactions_tab = ReactionsTab()
        self.rev_eng_modify_tab = ReverseEngineeringModifyTab()
        tabs.addTab(self.species_tab, "Species")
        tabs.addTab(self.reactions_tab, "Reactions")
        tabs.addTab(self.rev_eng_modify_tab, "Reverse Engineering")

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

        deterministic.triggered.connect(lambda _: DeterministicSimulationDialog())
        stochastic.triggered.connect(lambda _: StochasticSimulationDialog())

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

            self.species_tab.update_list()
            self.reactions_tab.update_list()


    # def _handler_remove_reactions_button(self):
    #     del GeneController.get_instance().network.reactions[self.reactions_panel.m_list.currentRow()]
    #     self._refresh_reactions_list(self.reactions_panel.m_list)
    #
    # def _handler_remove_species_button(self):
    #     species_id = self.species_panel.m_list.property("id" + str(self.species_panel.m_list.currentRow()))
    #     del GeneController.get_instance().network.species[species_id]
    #     self._refresh_species_list(self.species_panel.m_list)


app = QApplication([])
g = GeneWindow()
app.exec_()
