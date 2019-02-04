from models.network import Network


class GeneController:
    instance = None

    def __init__(self):
        self.network = Network()
        self.symbols = dict()

    @staticmethod
    def get_instance():
        if not GeneController.instance:
            GeneController.instance = GeneController()
        return GeneController.instance
