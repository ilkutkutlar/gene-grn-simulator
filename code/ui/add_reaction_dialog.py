from typing import List

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QVBoxLayout, QComboBox, QPushButton

from models.formulae import Formulae
from ui.gene_controller import GeneController
from models.reaction import Reaction
from ui.gui import validate_species


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
        self.transcribed_species = add_field_to_form("Transcribed Species")
        self.transcribed_species.setVisible(True)

        self.transcription_rate_field = add_number_field_to_form("Transcription Rate")
        self.transcription_rate_field.setVisible(True)

        self.kd_field = add_number_field_to_form("Kd")
        self.kd_field.setVisible(True)

        self.hill_coefficient_field = add_number_field_to_form("Hill coefficient")
        self.hill_coefficient_field.setVisible(True)

        # Translation
        self.translated_mrna = add_field_to_form("Translated mRNA Species")
        self.produced_protein = add_field_to_form("Produced Protein Species")
        self.translation_rate_field = add_number_field_to_form("Tranlation Rate")

        # mRNA and Protein decay
        self.decaying_species = add_number_field_to_form("Decaying Species")
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
        self.transcribed_species.setVisible(False)
        self.kd_field.setVisible(False)
        self.transcription_rate_field.setVisible(False)

        # Translation
        self.translation_rate_field.setVisible(False)
        self.translated_mrna.setVisible(False)
        self.produced_protein.setVisible(False)

        # Degradation
        self.decay_rate_field.setVisible(False)
        self.decaying_species.setVisible(False)

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
            if self.show_erro_message(message_text):
                return False
        else:
            return True

    def _handler_reaction_type_changed(self, index):
        self._hide_all_fields()

        if index == 0:
            self.transcription_rate_field.setVisible(True)
            self.kd_field.setVisible(True)
            self.hill_coefficient_field.setVisible(True)
            self.transcribed_species.setVisible(True)
        elif index == 1:
            self.translation_rate_field.setVisible(True)
            self.translated_mrna.setVisible(True)
            self.produced_protein.setVisible(True)
        elif index == 2 or index == 3:
            self.decay_rate_field.setVisible(True)
            self.decaying_species.setVisible(True)
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
            transcribed_species: str = self.transcribed_species.text().strip()

            if not self._error_check_species(transcribed_species):
                return

            r = Reaction([], [transcribed_species],
                         Formulae.transcription_rate(tr_rate, kd, hill_coeff, transcribed_species))
        elif index == 1:  # Translation
            tr_rate: float = to_float(0, self.translation_rate_field.text().strip())
            translated_mrna: str = self.translated_mrna.text().strip()
            produced_protein: str = self.produced_protein.text().strip()

            r = Reaction([translated_mrna], [produced_protein], Formulae.translation_rate(tr_rate, translated_mrna))
        elif index == 2 or index == 3:  # Degradation
            decay_rate: float = to_float(0, self.decay_rate_field.text().strip())
            decaying_species: str = self.decaying_species.text().strip()

            r = Reaction([decaying_species], [], Formulae.degradation_rate(decay_rate, decaying_species))
        else:  # index = 4, Custom reaction
            left_text = self.reactants_field.text().strip()
            right_text = self.products_field.text().strip()

            left: List[str] = [] if len(left_text) == 0 else left_text.split(",")
            right: List[str] = [] if len(right_text) == 0 else right_text.split(",")

            if not self._error_check_species(left, right):
                return

            eq = self.custom_equation_field.text().strip()
            r = Reaction(left, right, Formulae.custom_reaction_rate(eq))

        GeneController.get_instance().network.reactions.append(r)

        self.close()
