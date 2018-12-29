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
        self.combo.currentIndexChanged.connect(self._handler_reaction_type_changed)

    def _init_form_fields(self):
        self.fields = {}

        # TODO: Regulation!
        # Common
        self.fields["name"] = QLineEdit()
        self.fields["name"].setPlaceholderText("Name")

        # Transcription
        self.fields["transcription_rate"] = QLineEdit()
        self.fields["transcription_rate"].setPlaceholderText("Transcription Rate")

        self.fields["kd"] = QLineEdit()
        self.fields["kd"].setPlaceholderText("Kd")

        self.fields["hill_coefficient"] = QLineEdit()
        self.fields["hill_coefficient"].setPlaceholderText("Hill coefficient")

        # Translation
        self.fields["translation_rate"] = QLineEdit()
        self.fields["translation_rate"].setPlaceholderText("Tranlation Rate")

        # mRNA decay
        self.fields["mrna_decay_rate"] = QLineEdit()
        self.fields["mrna_decay_rate"].setPlaceholderText("Decay Rate")

        # Protein decay
        self.fields["protein_decay_rate"] = QLineEdit()
        self.fields["protein_decay_rate"].setPlaceholderText("Decay Rate")

        # Custom reaction
        self.fields["custom_equation"] = QLineEdit()
        self.fields["custom_equation"].setPlaceholderText("Equation")

        for field in self.fields:
            self.fields[field].setVisible(False)
            self.form.addWidget(self.fields[field])

        # Make the initial fields visible (i.e. the combo box has this reaction
        # selected when the dialog is first opened)
        self.fields["name"].setVisible(True)
        self.fields["transcription_rate"].setVisible(True)
        self.fields["kd"].setVisible(True)
        self.fields["hill_coefficient"].setVisible(True)

    def _handler_reaction_type_changed(self, index):
        for field in self.fields:
            self.fields[field].setVisible(False)

        self.fields["name"].setVisible(True)

        if index == 0:
            self.fields["transcription_rate"].setVisible(True)
            self.fields["kd"].setVisible(True)
            self.fields["hill_coefficient"].setVisible(True)
        elif index == 1:
            self.fields["translation_rate"].setVisible(True)
        elif index == 2:
            self.fields["mrna_decay_rate"].setVisible(True)
        elif index == 3:
            self.fields["protein_decay_rate"].setVisible(True)
        elif index == 4:
            self.fields["custom_equation"].setVisible(True)

    def _init_ui(self):
        self.layout = QVBoxLayout()  # The main layout
        self.form = QVBoxLayout()  # This holds the text fields for reaction parameters
        self.combo = QComboBox()  # Choose the type of reaction
        self.ok_button = QPushButton("OK", self)

        self._init_combo()
        self._init_form_fields()

        self.layout.addWidget(self.combo)
        self.layout.addLayout(self.form)
        self.layout.addWidget(self.ok_button)

        self.setLayout(self.layout)

        self.setWindowTitle("Add reaction")
        self.setWindowModality(Qt.WindowModal)
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
