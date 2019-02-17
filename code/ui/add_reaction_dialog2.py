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

        l1 = QLabel()
        l1.setText("Translation rate: ")
        self.translation_rate = QLineEdit()
        self.translation_fields.addRow(l1, self.translation_rate)

        l2 = QLabel()
        l2.setText("Translated mRNA: ")
        self.translated_mrna = QLineEdit()
        self.translation_fields.addRow(l2, self.translated_mrna)

        l3 = QLabel()
        l3.setText("Produced protein: ")
        self.produced_protein = QLineEdit()
        self.translation_fields.addRow(l3, self.produced_protein)

    def _init_degradation_fields(self):
        self.degradation_fields = QFormLayout()

        l1 = QLabel()
        l1.setText("Decay rate: ")
        self.decay_rate = QLineEdit()
        self.degradation_fields.addRow(l1, self.decay_rate)

        l2 = QLabel()
        l2.setText("Decaying species: ")
        self.decaying_species = QLineEdit()
        self.translation_fields.addRow(l2, self.decaying_species)

    def _init_custom_reaction_fields(self):
        self.custom_reaction_fields = QFormLayout()

        l1 = QLabel()
        l1.setText("Reactant: ")
        self.reactant = QLineEdit()
        self.custom_reaction_fields.addRow(l1, self.decay_rate)

        l2 = QLabel()
        l2.setText("Product: ")
        self.product = QLineEdit()
        self.custom_reaction_fields.addRow(l2, self.decaying_species)

        l3 = QLabel()
        l3.setText("Equation: ")
        self.equation = QLineEdit()
        self.custom_reaction_fields.addRow(l3, self.equation)

    def _init_transcription_reaction_fields(self):
        self.transcription_reaction_fields = QFormLayout()

        l1 = QLabel("Transcription rate")
        self.transcription_rate = QLineEdit()
        self.transcription_reaction_fields.addRow(l1, self.transcription_rate)

        l2 = QLabel("Hill coefficient: ")
        self.hill = QLineEdit()
        self.transcription_reaction_fields.addRow(l2, self.hill)

        l3 = QLabel("Kd: ")
        self.kd = QLineEdit()
        self.transcription_reaction_fields.addRow(l3, self.kd)

        l4 = QLabel("Transcribed species: ")
        self.transcribed_species = QLineEdit()
        self.transcription_reaction_fields.addRow(l4, self.transcribed_species)

        self.is_regulated = QCheckBox("Is Regulated")
        self.transcription_reaction_fields.addRow(self.is_regulated)

        self.activation_radio = QRadioButton("Activation")
        self.repression_radio = QRadioButton("Repression")
        self.activation_radio.setChecked(True)
        self.transcription_reaction_fields.addRow(self.activation_radio, self.repression_radio)

        self.regulator = QLineEdit()
        self.transcription_reaction_fields.addRow(self.regulator)
