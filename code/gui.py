import re
import sys
from typing import List

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, \
    QInputDialog, QDialog, QComboBox, QFormLayout, QLineEdit, QMainWindow, QAction, QMessageBox

from gene_controller import GeneController
from gillespie import GillespieSimulator
from models import TranslationReaction, TranscriptionReaction, MrnaDegradationReaction, ProteinDegradationReaction, \
    CustomReaction, SimulationSettings, Regulation, RegType


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

    def _validate_species(self, species):
        if species:
            for s in species:
                if s not in GeneController.get_instance().network.species:
                    return False
            return True
        else:
            return True

    def _error_check_species(self, left, right):
        error_message = QMessageBox()
        error_message.setIcon(QMessageBox.Warning)
        error_message.setWindowTitle("Error")
        error_message.setStandardButtons(QMessageBox.Ok)
        message_text = None

        if not self._validate_species(left):
            message_text = \
                "Reactants include some invalid species. " \
                "Please add the species before you use them."

        if not self._validate_species(right):
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

        if index == 0:
            GeneController.get_instance().network.reactions.append(
                TranscriptionReaction(self.transcription_rate_field.text(),
                                      self.kd_field.text(),
                                      self.hill_coefficient_field.text(), left=left, right=right))
        elif index == 1:
            GeneController.get_instance().network.reactions.append(
                TranslationReaction(self.translation_rate_field.text(), left=left, right=right)
            )
        elif index == 2:
            GeneController.get_instance().network.reactions.append(
                MrnaDegradationReaction(self.mrna_decay_rate_field.text(), left=left, right=right)
            )
        elif index == 3:
            GeneController.get_instance().network.reactions.append(
                ProteinDegradationReaction(self.protein_decay_rate_field.text(), left=left, right=right)
            )
        elif index == 4:
            GeneController.get_instance().network.reactions.append(
                CustomReaction(self.custom_equation_field.text(), left=left, right=right)
            )

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
        # self.remove_button.clicked.connect(remove_function)

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
            None)
        self.reactions_panel = AddRemoveListLayout(
            "Reactions",
            self._refresh_reactions_list,
            self._handler_add_reactions_button,
            None)
        self.regulations_panel = AddRemoveListLayout(
            "Regulations",
            self._refresh_regulations_list,
            self._handler_add_regulation_button,
            None)

        layout.addLayout(self.species_panel)
        layout.addLayout(self.reactions_panel)
        layout.addLayout(self.regulations_panel)
        self._init_menubar()

        central_widget = QWidget()
        central_widget.setLayout(layout)

        self.setCentralWidget(central_widget)
        self.setWindowTitle("Gene")
        self.show()

    def _handler_deterministic_clicked(self):
        print("Not implemented yet")

    def _handler_stochastic_clicked(self):
        (time, ok) = QInputDialog.getText(self, 'Simulation Settings', 'Simulation Time (s)')

        if ok:
            net = GeneController.get_instance().network
            labels = []
            for species in net.species:
                labels.append((species, species))

            s = SimulationSettings("Results", "Time", "Concentration", 0, time, labels)
            GillespieSimulator.visualise(GillespieSimulator.simulate(net, s), s)

    def _init_menubar(self):
        self.menubar = self.menuBar()
        file = self.menubar.addMenu("File")
        file.addAction("Open SBML file")
        simulate = self.menubar.addMenu("Simulate")

        deterministic = QAction("Deterministic (ODE)", self)
        stochastic = QAction("Stochastic (Gillespie Algorithm)", self)

        deterministic.triggered.connect(self._handler_deterministic_clicked)
        stochastic.triggered.connect(self._handler_stochastic_clicked)

        simulate.addAction(deterministic)
        simulate.addAction(stochastic)

    def _handler_add_reactions_button(self):
        dialog = AddReactionDialog()
        dialog.finished.connect(lambda: self._refresh_reactions_list(self.reactions_panel.m_list))
        dialog.exec_()

    def _handler_add_species_button(self):
        (species, ok) = QInputDialog.getText(self, 'Add new species', 'Species name:')

        if ok:
            GeneController.get_instance().network.species[species] = 0
            self._refresh_species_list(self.species_panel.m_list)

    def _handler_add_regulation_button(self):
        (reg, ok) = QInputDialog.getText(self, 'Add regulation', 'Format: X -> Y (Activator), X -| Y (Repressor)')

        if ok:
            reg_match = re.match(
                "(.*)(->|-\|)(.*)", reg)
            from_gene = reg_match.group(1)
            to_gene = reg_match.group(3)
            reg_type_str = reg_match.group(2)

            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Warning)
            error_message.setWindowTitle("Error")
            error_message.setStandardButtons(QMessageBox.Ok)
            message_text = "Regulation includes some invalid species. " \
                           "Please add the species before you use them."

            if not self._validate_species(from_gene) or not self._validate_species(to_gene):
                error_message.setText(message_text)

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

    @staticmethod
    def _refresh_reactions_list(m_list):
        m_list.clear()
        for reaction in GeneController.get_instance().network.reactions:
            m_list.addItem(reaction.__str__())

    @staticmethod
    def _refresh_species_list(m_list):
        m_list.clear()
        for species in GeneController.get_instance().network.species:
            m_list.addItem(species)

    @staticmethod
    def _refresh_regulations_list(m_list):
        m_list.clear()
        for regulation in GeneController.get_instance().network.regulations:
            m_list.addItem(regulation.__str__())


app = QApplication([])
g = GeneWindow()
app.exec_()
