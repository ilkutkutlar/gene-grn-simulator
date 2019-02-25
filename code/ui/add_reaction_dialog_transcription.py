from PyQt5.QtWidgets import QFormLayout, QRadioButton, QGroupBox, QLabel, QLineEdit, QWidget, QComboBox, QCheckBox, \
    QListWidget, QPushButton

from ui import common_widgets


class TranscriptionFields(QWidget):
    def __init__(self):
        super().__init__()
        fields = QFormLayout()

        self.reaction_name1 = QLineEdit()
        self.transcription_rate = QLineEdit()
        self.transcribed_species = common_widgets.make_species_combo()

        self._init_regulation()

        fields.addRow(QLabel("Reaction name"), self.reaction_name1)
        fields.addRow(QLabel("Transcription rate"), self.transcription_rate)
        fields.addRow(QLabel("Transcribed species: "), self.transcribed_species)

        fields.addRow(self.is_regulated)
        fields.addRow(self.hill_label, self.hill)
        fields.addRow(self.input_gate_label, self.input_gate)
        fields.addRow(self.regulations_label)
        fields.addRow(self.regulations_list)
        fields.addRow(self.add_regulation_box)

        self._is_regulated_state_changed()
        self.activation_radio.setChecked(True)
        self.setLayout(fields)

    def _add_regulation_clicked(self):
        pass

    def _refresh_regulations_list(self):
        self.regulations_list.clear()
        for x in self.regulations:
            self.regulations_list.addItem(x)

    def _is_regulated_state_changed(self):
        if self.is_regulated.isChecked():
            self.hill_label.setDisabled(False)
            self.hill.setDisabled(False)
            self.activation_radio.setDisabled(False)
            self.repression_radio.setDisabled(False)
            self.regulator.setDisabled(False)
            self.reg_label.setDisabled(False)
            self.k_label.setDisabled(False)
            self.k.setDisabled(False)
            self.regulations_label.setDisabled(False)
            self.input_gate_label.setDisabled(False)
            self.input_gate.setDisabled(False)
        else:
            self.hill_label.setDisabled(True)
            self.hill.setDisabled(True)
            self.activation_radio.setDisabled(True)
            self.repression_radio.setDisabled(True)
            self.regulator.setDisabled(True)
            self.reg_label.setDisabled(True)
            self.k_label.setDisabled(True)
            self.k.setDisabled(True)
            self.regulations_label.setDisabled(True)
            self.input_gate_label.setDisabled(True)
            self.input_gate.setDisabled(True)

    def _init_regulation(self):
        self.is_regulated = QCheckBox("Is Regulated")
        self.is_regulated.stateChanged.connect(self._is_regulated_state_changed)

        self.hill_label = QLabel("Hill coefficient: ")
        self.hill = QLineEdit()

        self.input_gate_label = QLabel("Input gate: ")
        self.input_gate = QComboBox()
        self.input_gate.addItem("AND")
        self.input_gate.addItem("OR")
        self.input_gate.addItem("SUM")

        self._init_add_regulation_box()

        self.regulations_label = QLabel("Regulations")
        self.regulations_list = QListWidget()
        self.regulations_list.setFixedHeight(100)

    def _init_add_regulation_box(self):
        self.add_regulation_box = QGroupBox()
        fields = QFormLayout()

        self.activation_radio = QRadioButton("Activation")
        self.repression_radio = QRadioButton("Repression")

        self.reg_label = QLabel("Regulating species:")
        self.regulator = common_widgets.make_species_combo()

        self.k_label = QLabel("K:")
        self.k = QLineEdit()

        self.add_button = QPushButton("Add regulation")
        self.add_button.clicked.connect(self._add_regulation_clicked)

        fields.addRow(self.activation_radio, self.repression_radio)
        fields.addRow(self.reg_label, self.regulator)
        fields.addRow(self.k_label, self.k)
        fields.addRow(self.add_button)
        self.add_regulation_box.setLayout(fields)
