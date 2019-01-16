import re
from typing import List

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, \
    QInputDialog, QDialog, QComboBox, QLineEdit, QMainWindow, QAction, QMessageBox, QFileDialog

from gene_controller import GeneController
from gillespie import GillespieSimulator
from models import TranslationReaction, TranscriptionReaction, MrnaDegradationReaction, ProteinDegradationReaction, \
    CustomReaction, SimulationSettings, Regulation, RegType, Network
from sbml_parser import SbmlParser


def validate_species(species):
    if species:
        for s in species:
            if s not in GeneController.get_instance().network.species:
                return False
        return True
    else:
        return True


class AddReactionDialog(QDialog):
    def __init__(self):
        super().__init__()
        self._init_ui()

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

        self.rp_info_field = QLabel()
        self.rp_info_field.setText("Reactants and products must be comma separated names of species")
        self.form.addWidget(self.rp_info_field)
        self.rp_info_field.setVisible(True)

        self.reactants_field = add_field_to_form("Reactants")
        self.reactants_field.setVisible(True)

        self.products_field = add_field_to_form("Products")
        self.products_field.setVisible(True)

        # Transcription
        self.transcription_rate_field = add_number_field_to_form("Transcription Rate")
        self.transcription_rate_field.setVisible(True)

        self.kd_field = add_number_field_to_form("Kd")
        self.kd_field.setVisible(True)

        self.hill_coefficient_field = add_number_field_to_form("Hill coefficient")
        self.hill_coefficient_field.setVisible(True)

        # Translation
        self.translation_rate_field = add_number_field_to_form("Tranlation Rate")

        # mRNA decay
        self.mrna_decay_rate_field = add_number_field_to_form("Decay Rate")

        # Protein decay
        self.protein_decay_rate_field = add_number_field_to_form("Decay Rate")

        # Custom reaction
        self.custom_equation_field = add_field_to_form("Equation")

    def _handler_reaction_type_changed(self, index):
        self.name_field.setVisible(True)
        self.rp_info_field.setVisible(True)
        self.reactants_field.setVisible(True)
        self.products_field.setVisible(True)

        self.transcription_rate_field.setVisible(False)
        self.kd_field.setVisible(False)
        self.hill_coefficient_field.setVisible(False)
        self.translation_rate_field.setVisible(False)
        self.mrna_decay_rate_field.setVisible(False)
        self.protein_decay_rate_field.setVisible(False)
        self.custom_equation_field.setVisible(False)

        if index == 0:
            self.transcription_rate_field.setVisible(True)
            self.kd_field.setVisible(True)
            self.hill_coefficient_field.setVisible(True)
        elif index == 1:
            self.translation_rate_field.setVisible(True)
        elif index == 2:
            self.mrna_decay_rate_field.setVisible(True)
        elif index == 3:
            self.protein_decay_rate_field.setVisible(True)
        elif index == 4:
            self.rp_info_field.setVisible(False)
            self.reactants_field.setVisible(False)
            self.products_field.setVisible(False)
            self.custom_equation_field.setVisible(True)

    def _error_check_species(self, left, right):
        error_message = QMessageBox()
        error_message.setIcon(QMessageBox.Warning)
        error_message.setWindowTitle("Error")
        error_message.setStandardButtons(QMessageBox.Ok)
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
            error_message.setText(message_text)
            button = error_message.exec_()
            if button == QMessageBox.Ok:
                return False
        else:
            return True

    def _handler_ok_button_clicked(self):
        index = self.combo.currentIndex()
        left_text = self.reactants_field.text().strip()
        right_text = self.products_field.text().strip()

        left: List[str] = [] if len(left_text) == 0 else left_text.split(",")
        right: List[str] = [] if len(right_text) == 0 else right_text.split(",")

        if not self._error_check_species(left, right):
            return

        def to_float(default, val) -> float:
            try:
                m: float = float(val)
            except ValueError:
                m: float = default
            return m

        if index == 0:
            tr_rate: float = to_float(
                0, self.transcription_rate_field.text().strip())

            kd: float = to_float(
                1, self.kd_field.text().strip())

            hill_coeff: float = to_float(
                1, self.hill_coefficient_field.text().strip())

            GeneController.get_instance().network.reactions.append(
                TranscriptionReaction(tr_rate, kd, hill_coeff, left=left, right=right))
        elif index == 1:
            tr_rate: float = to_float(0, self.translation_rate_field.text().strip())

            GeneController.get_instance().network.reactions.append(
                TranslationReaction(tr_rate, left=left, right=right))
        elif index == 2:
            decay_rate: float = to_float(0, self.mrna_decay_rate_field.text().strip())

            GeneController.get_instance().network.reactions.append(
                MrnaDegradationReaction(decay_rate, left=left, right=right))
        elif index == 3:
            decay_rate: float = to_float(0, self.protein_decay_rate_field.text().strip())

            GeneController.get_instance().network.reactions.append(
                ProteinDegradationReaction(decay_rate, left=left, right=right))
        elif index == 4:
            eq = self.custom_equation_field.text().strip()

            GeneController.get_instance().network.reactions.append(
                CustomReaction(eq, left=left, right=right))

        self.close()

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


class AddRemoveListLayout(QVBoxLayout):

    def __init__(self, label, refresh_function, add_function, remove_function):
        super().__init__()

        self.m_label = QLabel()
        self.m_label.setText(label)

        self.m_list = QListWidget()
        refresh_function(self.m_list)

        self.m_add_button = QPushButton("Add")
        self.m_remove_button = QPushButton("Remove")

        self.m_add_button.clicked.connect(add_function)
        self.m_remove_button.clicked.connect(remove_function)

        self.addWidget(self.m_label)
        self.addWidget(self.m_list)
        self.addWidget(self.m_add_button)
        self.addWidget(self.m_remove_button)
        self.addStretch()


class GeneWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout()

        self.species_panel = AddRemoveListLayout(
            "Species",
            self._refresh_species_list,
            self._handler_add_species_button,
            self._handler_remove_species_button)
        self.reactions_panel = AddRemoveListLayout(
            "Reactions",
            self._refresh_reactions_list,
            self._handler_add_reactions_button,
            self._handler_remove_reactions_button)
        self.regulations_panel = AddRemoveListLayout(
            "Regulations",
            self._refresh_regulations_list,
            self._handler_add_regulation_button,
            self._handler_remove_regulation_button)

        layout.addLayout(self.species_panel)
        layout.addLayout(self.reactions_panel)
        layout.addLayout(self.regulations_panel)
        self._init_menubar()

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
            self._refresh_regulations_list(self.regulations_panel.m_list)

    def _handler_deterministic_clicked(self):
        print("Not implemented yet")

    def _handler_stochastic_clicked(self):
        (time, ok) = QInputDialog.getText(self, 'Simulation Settings', 'Simulation Time (s)')

        if ok:
            end_time = int(time) if time else 0

            net = GeneController.get_instance().network
            labels = []
            for species in net.species:
                labels.append((species, species))

            s = SimulationSettings("Results", "Time", "Concentration", 0, end_time, labels)
            GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)

    def _handler_add_reactions_button(self):
        dialog = AddReactionDialog()
        dialog.finished.connect(lambda: self._refresh_reactions_list(self.reactions_panel.m_list))
        dialog.exec_()

    def _handler_remove_reactions_button(self):
        del GeneController.get_instance().network.reactions[self.reactions_panel.m_list.currentRow()]
        self._refresh_reactions_list(self.reactions_panel.m_list)

    def _handler_add_species_button(self):
        (species, ok) = QInputDialog.getText(self, 'Add new species', 'Species name:')

        if ok:
            GeneController.get_instance().network.species[species] = 0
            self._refresh_species_list(self.species_panel.m_list)

    def _handler_remove_species_button(self):
        species_id = self.species_panel.m_list.property("id" + str(self.species_panel.m_list.currentRow()))
        del GeneController.get_instance().network.species[species_id]
        self._refresh_species_list(self.species_panel.m_list)

    def _handler_add_regulation_button(self):
        (reg, ok) = QInputDialog.getText(self, 'Add regulation', 'Format: X -> Y (Activator), X -| Y (Repressor)')

        if ok:
            reg_match = re.match(
                "(.*)(->|-\|)(.*)", reg)
            from_gene = reg_match.group(1).strip()
            to_gene = reg_match.group(3).strip()
            reg_type_str = reg_match.group(2).strip()  # -> or -|

            if not validate_species(from_gene) or not validate_species(to_gene):
                error_message = QMessageBox()
                error_message.setIcon(QMessageBox.Warning)
                error_message.setWindowTitle("Error")
                error_message.setStandardButtons(QMessageBox.Ok)
                error_message.setText("Regulation includes some invalid species. "
                                      "Please add the species before you use them.")

                button = error_message.exec_()
                if button == QMessageBox.Ok:
                    return

            reg_type = None
            if reg_type_str == "->":
                reg_type = RegType.ACTIVATION
            elif reg_type_str == "-|":
                reg_type = RegType.REPRESSION
            else:
                print("ERROR!")

            GeneController.get_instance().network.regulations.append(Regulation(from_gene, to_gene, reg_type))
            self._refresh_regulations_list(self.regulations_panel.m_list)

    def _handler_remove_regulation_button(self):
        del GeneController.get_instance().network.regulations[self.regulations_panel.m_list.currentRow()]
        self._refresh_reactions_list(self.regulations_panel.m_list)

    @staticmethod
    def _refresh_reactions_list(m_list):
        m_list.clear()
        for reaction in GeneController.get_instance().network.reactions:
            m_list.addItem(reaction.__str__())
        print(GeneController.get_instance().network.reactions)

    @staticmethod
    def _refresh_species_list(m_list):
        m_list.clear()
        i: int = 0
        for species in GeneController.get_instance().network.species:
            m_list.addItem(species)
            m_list.setProperty("id" + str(i), species)
            i += 1
        print(GeneController.get_instance().network.species)

    @staticmethod
    def _refresh_regulations_list(m_list):
        m_list.clear()
        for regulation in GeneController.get_instance().network.regulations:
            m_list.addItem(regulation.__str__())
        print(GeneController.get_instance().network.regulations)


app = QApplication([])
g = GeneWindow()
app.exec_()
