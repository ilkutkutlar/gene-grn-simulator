from PyQt5.QtWidgets import QListWidget, QLabel, QGridLayout, QPushButton, \
    QVBoxLayout, QDialog, QFormLayout, QLineEdit, QWidget, QCheckBox, QRadioButton, QComboBox

from models.formulae import TranscriptionFormula, TranslationFormula, DegradationFormula, CustomFormula
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation
from ui.gene_controller import GeneController


class AddReactionDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()
        self.grid = QGridLayout()
        self.properties = QWidget()
        self.properties_layout = QVBoxLayout()
        self.properties.setLayout(self.properties_layout)

        self._init_transcription_fields()
        self._init_translation_fields()
        self._init_degradation_fields()
        self._init_custom_reaction_fields()

        self._init_reaction_types()
        self.grid.addWidget(self.properties, 0, 1)

        self.main_layout.addLayout(self.grid)
        self.setLayout(self.main_layout)
        self._init_ok_button()
        self.setWindowTitle("Add Reaction")

    def _make_species_combo(self):
        combo = QComboBox()
        for x in GeneController.get_instance().get_species():
            combo.addItem(x)
        return combo

    def _reaction_type_changed_handler(self):
        index = self.reaction_types_list.currentRow()

        self.transcription_fields.hide()
        self.translation_fields.hide()
        self.degradation_fields.hide()
        self.custom_fields.hide()

        if index == 0:  # Transcription
            self.transcription_fields.show()
        elif index == 1:  # Translation
            self.translation_fields.show()
        elif index == 2:  # Degradation
            self.degradation_fields.show()
        else:  # Custom
            self.custom_fields.show()

    def _ok_button_clicked_handler(self):
        index = self.reaction_types_list.currentRow()

        if index == 0:  # Transcription
            name = str(self.reaction_name1.text())
            rate = float(self.transcription_rate.text())
            n = float(self.hill.text())
            kd = float(self.kd.text())
            species = str(self.transcribed_species.currentText())

            regs = []
            if self.is_regulated.isChecked():
                reg_type = RegType.ACTIVATION if self.activation_radio.isChecked() else RegType.REPRESSION
                regs.append(Regulation(self.regulator.currentText(), species, reg_type))
            formula = TranscriptionFormula(rate, n, kd, species, regs)
            GeneController.get_instance().add_reaction(Reaction(name, [], [species], formula))

        elif index == 1:  # Translation
            name = str(self.reaction_name2.text())
            rate = float(self.translation_rate.text())
            translated_mrna = str(self.translated_mrna.currentText())
            produced_protein = str(self.produced_protein.currentText())
            formula = TranslationFormula(rate, translated_mrna)
            GeneController.get_instance().add_reaction(Reaction(name, [translated_mrna], [produced_protein], formula))

        elif index == 2:  # Degradation
            name = str(self.reaction_name3.text())
            rate = float(self.decay_rate.text())
            decaying_species = str(self.decaying_species.currentText())
            formula = DegradationFormula(rate, decaying_species)
            GeneController.get_instance().add_reaction(Reaction(name, [decaying_species], [], formula))
        else:  # Custom
            name = str(self.reaction_name4.text())
            reactant = str(self.reactant.currentText())
            product = str(self.product.currentText())
            equation = str(self.equation.text())
            formula = CustomFormula(equation, {}, {})
            GeneController.get_instance().add_reaction(Reaction(name, [reactant], [product], formula))

        self.close()

    def _init_reaction_types(self):
        self.reaction_types_list = QListWidget()

        self.reaction_types_list.setMinimumWidth(200)
        self.reaction_types_list.setMaximumWidth(200)

        self.reaction_types_list.addItem("Transcription Reaction")
        self.reaction_types_list.addItem("Translation Reaction")
        self.reaction_types_list.addItem("Degradation Reaction")
        self.reaction_types_list.addItem("Custom Reaction")

        self.grid.addWidget(self.reaction_types_list, 0, 0)
        self.reaction_types_list.itemClicked.connect(self._reaction_type_changed_handler)

    def _init_ok_button(self):
        self.ok_button = QPushButton("Add reaction")
        self.ok_button.clicked.connect(self._ok_button_clicked_handler)
        self.main_layout.addWidget(self.ok_button)

    def _init_transcription_fields(self):
        fields = QFormLayout()

        self.reaction_name1 = QLineEdit()
        fields.addRow(QLabel("Reaction name"), self.reaction_name1)
        self.transcription_rate = QLineEdit()
        fields.addRow(QLabel("Transcription rate"), self.transcription_rate)
        self.hill = QLineEdit()
        fields.addRow(QLabel("Hill coefficient: "), self.hill)
        self.kd = QLineEdit()
        fields.addRow(QLabel("Kd: "), self.kd)
        self.transcribed_species = self._make_species_combo()
        fields.addRow(QLabel("Transcribed species: "), self.transcribed_species)

        def is_regulated_state_changed():
            if self.is_regulated.isChecked():
                self.activation_radio.setDisabled(False)
                self.repression_radio.setDisabled(False)
                self.regulator.setDisabled(False)
                self.reg_label.setDisabled(False)
            else:
                self.activation_radio.setDisabled(True)
                self.repression_radio.setDisabled(True)
                self.regulator.setDisabled(True)
                self.reg_label.setDisabled(True)

        self.is_regulated = QCheckBox("Is Regulated")
        self.is_regulated.stateChanged.connect(is_regulated_state_changed)

        self.activation_radio = QRadioButton("Activation")
        self.activation_radio.setDisabled(True)
        self.activation_radio.setChecked(True)

        self.repression_radio = QRadioButton("Repression")
        self.repression_radio.setDisabled(True)

        self.reg_label = QLabel("Regulating species:")
        self.regulator = self._make_species_combo()
        self.regulator.setDisabled(True)
        self.reg_label.setDisabled(True)

        fields.addRow(self.is_regulated)
        fields.addRow(self.activation_radio, self.repression_radio)
        fields.addRow(self.reg_label, self.regulator)

        self.transcription_fields = QWidget()
        self.transcription_fields.setLayout(fields)
        self.properties_layout.addWidget(self.transcription_fields)

    def _init_translation_fields(self):
        fields = QFormLayout()

        self.reaction_name2 = QLineEdit()
        fields.addRow(QLabel("Reaction name"), self.reaction_name2)
        self.translation_rate = QLineEdit()
        fields.addRow(QLabel("Translation rate: "), self.translation_rate)
        self.translated_mrna = self._make_species_combo()
        fields.addRow(QLabel("Translated mRNA: "), self.translated_mrna)
        self.produced_protein = self._make_species_combo()
        fields.addRow(QLabel("Produced protein: "), self.produced_protein)

        self.translation_fields = QWidget()
        self.translation_fields.setLayout(fields)
        self.properties_layout.addWidget(self.translation_fields)
        self.translation_fields.hide()

    def _init_degradation_fields(self):
        fields = QFormLayout()

        self.reaction_name3 = QLineEdit()
        fields.addRow(QLabel("Reaction name"), self.reaction_name3)
        self.decay_rate = QLineEdit()
        fields.addRow(QLabel("Decay rate: "), self.decay_rate)
        self.decaying_species = self._make_species_combo()
        fields.addRow(QLabel("Decaying species: "), self.decaying_species)

        self.degradation_fields = QWidget()
        self.degradation_fields.setLayout(fields)
        self.properties_layout.addWidget(self.degradation_fields)
        self.degradation_fields.hide()

    def _init_custom_reaction_fields(self):
        fields = QFormLayout()

        self.reaction_name4 = QLineEdit()
        fields.addRow(QLabel("Reaction name"), self.reaction_name4)
        self.equation = QLineEdit()
        fields.addRow(QLabel("Equation: "), self.equation)
        self.reactant = self._make_species_combo()
        fields.addRow(QLabel("Reactant: "), self.reactant)
        self.product = self._make_species_combo()
        fields.addRow(QLabel("Product: "), self.product)

        self.custom_fields = QWidget()
        self.custom_fields.setLayout(fields)
        self.properties_layout.addWidget(self.custom_fields)
        self.custom_fields.hide()
