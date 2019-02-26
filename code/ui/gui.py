from PyQt5.QtWidgets import QTabWidget
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QWidget, QApplication, QMainWindow, QAction, QMessageBox

from ui.simulation.deterministic_simulation_dialog import DeterministicSimulationDialog
from ui.open_sbml_dialog import OpenSbmlDialog
from ui.reactions.reactions_tab import ReactionsTab
from ui.reverse_engineering.reverse_engineering_modify_tab import ReverseEngineeringModifyTab
from ui.species.species_tab import SpeciesTab
from ui.simulation.stochastic_simulation_dialog import StochasticSimulationDialog


def show_error_message(message):
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
        open_file.triggered.connect(lambda _: OpenSbmlDialog(self))

        file.addAction(open_file)

        # Simulate menu
        simulate = self.menubar.addMenu("Simulate")

        deterministic = QAction("Deterministic (ODE)", self)
        deterministic.triggered.connect(lambda _: DeterministicSimulationDialog())
        simulate.addAction(deterministic)

        stochastic = QAction("Stochastic (Gillespie Algorithm)", self)
        stochastic.triggered.connect(lambda _: StochasticSimulationDialog())
        simulate.addAction(stochastic)


app = QApplication([])
g = GeneWindow()
app.exec_()
