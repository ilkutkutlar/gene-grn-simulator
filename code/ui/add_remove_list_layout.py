from PyQt5.QtWidgets import QVBoxLayout, QLabel, QListWidget, QPushButton


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
