from typing import List

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QComboBox, QTabWidget
from PyQt5.QtWidgets import QDialog, QLineEdit
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QListWidget, QPushButton
from PyQt5.QtWidgets import QWidget, QApplication, QHBoxLayout, QMainWindow, QAction, QMessageBox, \
    QFileDialog

from input_output.sbml_parser import SbmlParser
from models.formulae import TranscriptionFormula, TranslationFormula, DegradationFormula, CustomFormula
from simulation.ode_simulator import OdeSimulator
from ui.gene_controller import GeneController
from simulation.gillespie_simulator import GillespieSimulator
from models.network import Network
from models.reaction import Reaction
from models.simulation_settings import SimulationSettings
from ui.multiple_input_dialog import QMultipleInputDialog
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


class AddRemoveListLayout(QVBoxLayout):

    def __init__(self, label, refresh_function, add_function, remove_function):
        super().__init__()

        self.m_label = QLabel()
        self.m_label.setText(label)

        self.m_list = QListWidget()

        self.m_add_button = QPushButton("Add")
        self.m_remove_button = QPushButton("Remove")

        self.m_add_button.clicked.connect(add_function)
        self.m_remove_button.clicked.connect(remove_function)

        self.addWidget(self.m_label)
        self.addWidget(self.m_list)
        self.addStretch()
        self.addWidget(self.m_add_button)
        self.addWidget(self.m_remove_button)

        refresh_function(self.m_list)


class AddReactionDialog(QDialog):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        self.layout = QVBoxLayout()  # The main layout
        self.form = QVBoxLayout()  # This holds the text fields for reaction parameters
        self.combo = QComboBox()  # Choose the type of reaction
        self.ok_button = QPushButton("OK", self)
        self.ok_button.clicked.connect(self._handler_ok_button_clicked)

        self._init_combo()
        self._init_form_fields()

        self.layout.addWidget(self.combo)
        self.layout.addLayout(self.form)
        self.layout.addWidget(self.ok_button)

        self.setLayout(self.layout)

        self.setWindowTitle("Add reaction")
        self.setWindowModality(Qt.WindowModal)

    def _init_combo(self):
        self.combo.addItem("Transcription Reaction")
        self.combo.addItem("Translation Reaction")
        self.combo.addItem("mRNA Degradation Reaction")
        self.combo.addItem("Protein Degradation Reaction")
        self.combo.addItem("Custom Reaction")
        self.combo.currentIndexChanged.connect(self._handler_reaction_type_changed)

    def _init_form_fields(self):
        validator = QDoubleValidator()
        # Standard notation disallows constants such as "e"
        validator.setNotation(QDoubleValidator.StandardNotation)

        def add_field_to_form(placeholder):
            field = QLineEdit()
            field.setPlaceholderText(placeholder)
            field.setVisible(False)
            self.form.addWidget(field)
            return field

        def add_number_field_to_form(placeholder):
            field = add_field_to_form(placeholder)
            field.setValidator(validator)
            return field

        # Common
        # Make the initial fields visible (i.e. the combo box has this reaction
        # selected when the dialog is first opened)
        self.name_field = add_field_to_form("Name")
        self.name_field.setVisible(True)

        # Transcription
        self.transcribed_species_field = add_field_to_form("Transcribed Species")
        self.transcribed_species_field.setVisible(True)

        self.transcription_rate_field = add_number_field_to_form("Transcription Rate")
        self.transcription_rate_field.setVisible(True)

        self.kd_field = add_number_field_to_form("Kd")
        self.kd_field.setVisible(True)

        self.hill_coefficient_field = add_number_field_to_form("Hill coefficient")
        self.hill_coefficient_field.setVisible(True)

        # Translation
        self.translated_mrna_field = add_field_to_form("Translated mRNA Species")
        self.produced_protein_field = add_field_to_form("Produced Protein Species")
        self.translation_rate_field = add_number_field_to_form("Tranlation Rate")

        # mRNA and Protein decay
        self.decaying_species_field = add_number_field_to_form("Decaying Species")
        self.decay_rate_field = add_number_field_to_form("Decay Rate")

        # Custom reaction
        self.rp_info_field = QLabel()
        self.rp_info_field.setText("Reactants and products must be comma separated names of species")
        self.form.addWidget(self.rp_info_field)
        self.rp_info_field.setVisible(False)

        self.reactants_field = add_field_to_form("Reactants")
        self.products_field = add_field_to_form("Products")
        self.custom_equation_field = add_field_to_form("Equation")

    def _hide_all_fields(self):
        # Common
        self.name_field.setVisible(False)

        # Transcription
        self.hill_coefficient_field.setVisible(False)
        self.transcribed_species_field.setVisible(False)
        self.kd_field.setVisible(False)
        self.transcription_rate_field.setVisible(False)

        # Translation
        self.translation_rate_field.setVisible(False)
        self.translated_mrna_field.setVisible(False)
        self.produced_protein_field.setVisible(False)

        # Degradation
        self.decay_rate_field.setVisible(False)
        self.decaying_species_field.setVisible(False)

        # Custom reactions
        self.rp_info_field.setVisible(False)
        self.reactants_field.setVisible(False)
        self.products_field.setVisible(False)
        self.custom_equation_field.setVisible(False)

    def _error_check_species(self, left, right):
        message_text = None

        if not validate_species(left):
            message_text = \
                "Reactants include some invalid species. " \
                "Please add the species before you use them."

        if not validate_species(right):
            message_text = \
                "Products include some invalid species. " \
                "Please add the species before you use them."

        if message_text:
            if show_error_message(message_text):
                return False
        else:
            return True

    def _handler_reaction_type_changed(self, index):
        self._hide_all_fields()
        self.name_field.setVisible(True)

        if index == 0:
            self.transcription_rate_field.setVisible(True)
            self.kd_field.setVisible(True)
            self.hill_coefficient_field.setVisible(True)
            self.transcribed_species_field.setVisible(True)
        elif index == 1:
            self.translation_rate_field.setVisible(True)
            self.translated_mrna_field.setVisible(True)
            self.produced_protein_field.setVisible(True)
        elif index == 2 or index == 3:
            self.decay_rate_field.setVisible(True)
            self.decaying_species_field.setVisible(True)
        elif index == 4:
            self.rp_info_field.setVisible(True)
            self.reactants_field.setVisible(True)
            self.products_field.setVisible(True)
            self.custom_equation_field.setVisible(True)

    def _handler_ok_button_clicked(self):
        index = self.combo.currentIndex()

        def to_float(default, val) -> float:
            try:
                m: float = float(val)
            except ValueError:
                m: float = default
            return m

        if index == 0:  # Transcription
            tr_rate: float = to_float(0, self.transcription_rate_field.text().strip())
            kd: float = to_float(1, self.kd_field.text().strip())
            hill_coeff: float = to_float(1, self.hill_coefficient_field.text().strip())
            transcribed_species: str = self.transcribed_species_field.text().strip()

            if not self._error_check_species([], transcribed_species):
                return

            r = Reaction([], [transcribed_species],
                         TranscriptionFormula(tr_rate, kd, hill_coeff, transcribed_species))
        elif index == 1:  # Translation
            tr_rate: float = to_float(0, self.translation_rate_field.text().strip())
            translated_mrna: str = self.translated_mrna_field.text().strip()
            produced_protein: str = self.produced_protein_field.text().strip()

            if not self._error_check_species(translated_mrna, produced_protein):
                return

            r = Reaction([translated_mrna], [produced_protein], TranslationFormula(tr_rate, translated_mrna))
        elif index == 2 or index == 3:  # Degradation
            decay_rate: float = to_float(0, self.decay_rate_field.text().strip())
            decaying_species: str = self.decaying_species_field.text().strip()

            if not self._error_check_species(decaying_species, []):
                return

            r = Reaction([decaying_species], [], DegradationFormula(decay_rate, decaying_species))
        else:  # index = 4, Custom reaction
            left_text = self.reactants_field.text().strip()
            right_text = self.products_field.text().strip()

            left: List[str] = [] if len(left_text) == 0 else left_text.split(",")
            right: List[str] = [] if len(right_text) == 0 else right_text.split(",")

            if not self._error_check_species(left, right):
                return

            eq = self.custom_equation_field.text().strip()
            r = Reaction(left, right, CustomFormula(eq))

        GeneController.get_instance().network.reactions.append(r)

        self.close()


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

        info = QLabel()
        info.setText("Species are comma separated")

        time_field = QLineEdit()
        time_field.setPlaceholderText("Simulation time")

        sampling_field = QLineEdit()
        sampling_field.setPlaceholderText("Sampling rate")

        species_field = QLineEdit()
        species_field.setPlaceholderText("Which species to show")

        ok_button = QPushButton()
        ok_button.setText("Ok")

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

        info = QLabel()
        info.setText("Species are comma separated")

        time_field = QLineEdit()
        time_field.setPlaceholderText("Simulation time")

        species_field = QLineEdit()
        species_field.setPlaceholderText("Which species to show")

        ok_button = QPushButton()
        ok_button.setText("Ok")

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

    def _handler_add_reactions_button(self):
        dialog = AddReactionDialog()
        dialog.finished.connect(lambda: self._refresh_reactions_list(self.reactions_panel.m_list))
        dialog.exec_()

    def _handler_remove_reactions_button(self):
        del GeneController.get_instance().network.reactions[self.reactions_panel.m_list.currentRow()]
        self._refresh_reactions_list(self.reactions_panel.m_list)

    def _handler_add_species_button(self):
        dia = QDialog()
        main = QVBoxLayout()

        name = QLineEdit()
        name.setPlaceholderText("Species name")

        validator = QDoubleValidator()
        # Standard notation disallows constants such as "e"
        validator.setNotation(QDoubleValidator.StandardNotation)
        init_con = QLineEdit()
        init_con.setPlaceholderText("Initial amount")
        init_con.setValidator(validator)

        ok_button = QPushButton()
        ok_button.setText("Ok")

        def dialog_ok_button_clicked_handler():
            name_text = name.text().strip()
            init_con_text = init_con.text().strip()
            init_con_value = float(init_con_text) if init_con_text else float(0)

            if name_text == "":
                if show_error_message("You can't have a blank species name."):
                    return
            else:
                GeneController.get_instance().network.species[name_text] = init_con_value
                self._refresh_species_list(self.species_panel.m_list)
            dia.close()

        ok_button.clicked.connect(dialog_ok_button_clicked_handler)

        main.addWidget(name)
        main.addWidget(init_con)
        main.addWidget(ok_button)
        dia.setFixedHeight(120)
        dia.setMinimumWidth(200)
        dia.setWindowTitle("Add new species")
        dia.setLayout(main)
        dia.exec_()

    def _handler_remove_species_button(self):
        species_id = self.species_panel.m_list.property("id" + str(self.species_panel.m_list.currentRow()))
        del GeneController.get_instance().network.species[species_id]
        self._refresh_species_list(self.species_panel.m_list)

    @staticmethod
    def _refresh_reactions_list(m_list):
        m_list.clear()
        for reaction in GeneController.get_instance().network.reactions:
            m_list.addItem(reaction.__str__())

    @staticmethod
    def _refresh_species_list(m_list):
        m_list.clear()
        i: int = 0
        for species in GeneController.get_instance().network.species:
            m_list.addItem(species + ": " + str(GeneController.get_instance().network.species[species]))
            m_list.setProperty("id" + str(i), species)
            i += 1


app = QApplication([])
g = GeneWindow()
app.exec_()
