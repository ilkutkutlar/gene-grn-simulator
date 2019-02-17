from PyQt5.QtWidgets import QListWidget, QLabel, QGridLayout, QPushButton, \
    QVBoxLayout, QDialog, QFormLayout, QLineEdit, QWidget, QCheckBox, QRadioButton


class AddReactionDialog2(QDialog):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()
        self.grid = QGridLayout()
        self.reaction_properties_base = QWidget()

        self._init_transcription_reaction_fields()

        self._init_reaction_types_list()
        self.grid.addWidget(self.reaction_properties_base, 0, 1)
        self.reaction_properties_base.setLayout(self.transcription_reaction_fields)
        self._init_ok_button()

        self.main_layout.addLayout(self.grid)
        self.setLayout(self.main_layout)

    def _init_ok_button(self):
        self.ok_button = QPushButton("Add reaction")
        # self.add_button.clicked.connect(self._add_species_click_handler)
        self.main_layout.addWidget(self.ok_button)

    def _init_reaction_types_list(self):
        self.reaction_types_list = QListWidget()

        self.reaction_types_list.setMinimumWidth(200)
        self.reaction_types_list.setMaximumWidth(200)

        self.reaction_types_list.addItem("Transcription Reaction")
        self.reaction_types_list.addItem("Translation Reaction")
        self.reaction_types_list.addItem("Degradation Reaction")
        self.reaction_types_list.addItem("Custom Reaction")

        self.grid.addWidget(self.reaction_types_list, 0, 0)

    def _init_translation_fields(self):
        self.translation_fields = QFormLayout()

        self.translation_rate = QLineEdit()
        self.translation_fields.addRow(
            QLabel("Translation rate: "),
            self.translation_rate)

        self.translated_mrna = QLineEdit()
        self.translation_fields.addRow(
            QLabel("Translated mRNA: "),
            self.translated_mrna)

        self.produced_protein = QLineEdit()
        self.translation_fields.addRow(
            QLabel("Produced protein: "),
            self.produced_protein)

    def _init_degradation_fields(self):
        self.degradation_fields = QFormLayout()

        self.decay_rate = QLineEdit()
        self.degradation_fields.addRow(
            QLabel("Decay rate: "),
            self.decay_rate)

        self.decaying_species = QLineEdit()
        self.translation_fields.addRow(
            QLabel("Decaying species: "),
            self.decaying_species)

    def _init_custom_reaction_fields(self):
        self.custom_reaction_fields = QFormLayout()

        self.reactant = QLineEdit()
        self.custom_reaction_fields.addRow(
            QLabel("Reactant: "),
            self.decay_rate)

        self.product = QLineEdit()
        self.custom_reaction_fields.addRow(
            QLabel("Product: "),
            self.decaying_species)

        self.equation = QLineEdit()
        self.custom_reaction_fields.addRow(
            QLabel("Equation: "),
            self.equation)

    def _init_transcription_reaction_fields(self):
        self.transcription_reaction_fields = QFormLayout()

        self.transcription_rate = QLineEdit()
        self.transcription_reaction_fields.addRow(
            QLabel("Transcription rate"),
            self.transcription_rate)

        self.hill = QLineEdit()
        self.transcription_reaction_fields.addRow(
            QLabel("Hill coefficient: "),
            self.hill)

        self.kd = QLineEdit()
        self.transcription_reaction_fields.addRow(
            QLabel("Kd: "),
            self.kd)

        self.transcribed_species = QLineEdit()
        self.transcription_reaction_fields.addRow(
            QLabel("Transcribed species: "),
            self.transcribed_species)

        def is_regulated_state_changed():
            if self.is_regulated.isChecked():
                self.activation_radio.setDisabled(False)
                self.repression_radio.setDisabled(False)
                self.regulator.setDisabled(False)
            else:
                self.activation_radio.setDisabled(True)
                self.repression_radio.setDisabled(True)
                self.regulator.setDisabled(True)

        self.is_regulated = QCheckBox("Is Regulated")
        self.is_regulated.stateChanged.connect(is_regulated_state_changed)
        self.transcription_reaction_fields.addRow(self.is_regulated)

        self.activation_radio = QRadioButton("Activation")
        self.activation_radio.setDisabled(True)
        self.activation_radio.setChecked(True)
        self.repression_radio = QRadioButton("Repression")
        self.repression_radio.setDisabled(True)
        self.transcription_reaction_fields.addRow(self.activation_radio, self.repression_radio)

        self.regulator = QLineEdit()
        self.regulator.setDisabled(True)
        self.transcription_reaction_fields.addRow(self.regulator)
