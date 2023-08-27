from tasks import TaskHandler
from model.rnadata import RNAData
import utils


class UIConnection:

	def __init__(self) -> None:
		self.functions = {}
		self.ui_elemets = {}

	def add_ui_element(self, name, element) -> None:
		self.ui_elemets[name] = element

	def add_function(self, name, function) -> None:
		self.functions[name] = function

	def get_element(self, ui_element_name):
		return self.ui_elemets.get(ui_element_name)

	def call(self, function_name, *args, **kwargs) -> None:
		function = self.functions.get(function_name)
		if function:
			function(self, *args, **kwargs)


class Controller:

	class ModelInterface:

		def __init__(self, controller):
			self.controller = controller

		def __setattr__(self, key: str, value):
			if hasattr(self, 'controller'):
				print(f'setting {key}={value}')
				setattr(self.controller.current_task.model, key, value)
			else:
				super().__setattr__(key, value)

	def __init__(self, view) -> None:
		self.view = view
		self.current_task = None
		self.model_interface = self.ModelInterface(self)

	def query_model(self, query):
		if self.current_task:
			return getattr(self.current_task.model, query, None)
		return None

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		rna_data = RNAData(seq_record)
		task = TaskHandler(self, rna_data)
		self.current_task = task
		self.view.open_page('metainf', task_id=self.current_task.id)

	def revert_model_changes(self) -> None:
		self.current_task.reload_model()

	def start_task(self, options) -> None:
		self.current_task.start(options)
		self.view.open_page('results', task_id=self.current_task.id)

	def cancel_task(self) -> None:
		if self.current_task:
			self.current_task.cancel()
			self.current_task = None

	def connect_ui(self) -> UIConnection:
		ui_connection = UIConnection()
		self.current_task.connect_ui(ui_connection)
		return ui_connection
