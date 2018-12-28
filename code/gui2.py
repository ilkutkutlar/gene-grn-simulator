import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QApplication, QLabel, QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, \
    QInputDialog, QDialog, QComboBox, QFormLayout, QLineEdit


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

        self.combo.currentIndexChanged.connect(self._reaction_type_changed)
        self.combo.setCurrentIndex(0)

    def _init_fields(self):
        self.fields = {}

        # TODO: Regulation!
        # Common
        self.fields["name"] = QLineEdit()

        # Transcription
        self.edit_transcription_rate = QLineEdit()
        self.edit_kd = QLineEdit()
        self.edit_hill_coeff = QLineEdit()

        # Translation
        self.edit_translation_rate = QLineEdit()

        # mRNA decay
        self.edit_mrna_decay_rate = QLineEdit()

        # Protein decay
        self.edit_protein_decay_rate = QLineEdit()

        # Custom reaction
        self.edit_equation = QLineEdit()

    def _clear_form_fields(self):


def _init_ui(self):
        self.layout = QVBoxLayout()
        self.form = QFormLayout()
        self.combo = QComboBox()
        self.ok_button = QPushButton("OK", self)

        self._init_combo()
        self._init_fields()

        form.addRow("Name: ", name)
        self.form.addRow("Transcription Rate: ", trans_rate)
        self.form.addRow("Kd: ", kd)
        self.form.addRow("Hill coefficient: ", hill_coeff)
        form.addRow("Tranlation Rate: ", trans_rate)
        form.addRow("Decay Rate: ", decay_rate)
        form.addRow("Decay Rate: ", decay_rate)
        form.addRow("Equation: ", equation)

        self.layout.addWidget(self.combo)
        self.layout.addLayout(self.form)
        self.layout.addWidget(self.ok_button)

        self.setLayout(self.layout)
        self.setWindowTitle("Add reaction")
        self.setWindowModality(Qt.ApplicationModal)
        self.exec_()


class GeneWindow(QWidget):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout()
        layout.addLayout(self._species())
        layout.addLayout(self._reactions())

        self.setLayout(layout)
        self.setWindowTitle("Gene")
        self.show()

    def _add_reactions(self):
        dialog = AddReactionDialog()

    def _reactions(self):
        reactions_label = QLabel(self)
        reactions_label.setText("Reactions")
        reactions_list = QListWidget()
        reactions_list.addItem("Test")
        add_button = QPushButton("Add")
        remove_button = QPushButton("Remove")
        add_button.clicked.connect(self._add_reactions)

        layout = QVBoxLayout()
        layout.addWidget(reactions_label)
        layout.addWidget(reactions_list)
        layout.addWidget(add_button)
        layout.addWidget(remove_button)
        layout.addStretch()

        return layout

    def _add_species(self):
        (species, ok) = QInputDialog.getText(self, 'Add new species', 'Species name:')

        if ok:
            print(str(species))

    def _species(self):
        species_label = QLabel(self)
        species_label.setText("Species")
        species_list = QListWidget()
        species_list.addItem("Test")
        add_button = QPushButton("Add")
        remove_button = QPushButton("Remove")
        add_button.clicked.connect(self._add_species)

        layout = QVBoxLayout()
        layout.addWidget(species_label)
        layout.addWidget(species_list)
        layout.addWidget(add_button)
        layout.addWidget(remove_button)
        layout.addStretch()
        return layout


app = QApplication([])
g = GeneWindow()
app.exec_()
