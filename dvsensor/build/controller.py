from operator import attrgetter
from tasks import TaskHandler
from model.taskdata import SingleSeqAnalysisData
from model.rnadata import RNAData
import utils
import asyncio
import threading


class UIConnection:

	def __init__(self) -> None:
		self.functions = {}
		self.ui_elements = {}
		self.page_render_complete = threading.Event()

	def add_ui_element(self, name, element) -> None:
		self.ui_elements[name] = element

	def add_function(self, name, function) -> None:
		self.functions[name] = function

	def get_element(self, ui_element_name):
		return self.ui_elements.get(ui_element_name, None)

	def call(self, function_name, *args, **kwargs) -> None:
		function = self.functions.get(function_name)
		if function:
			function(self, *args, **kwargs)


class Controller:

	def __init__(self, view) -> None:
		self.view = view
		self.current_task = None

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		task_data = SingleSeqAnalysisData(rna_data=RNAData(seq_record))
		task = TaskHandler(self, task_data)
		self.current_task = task
		self.view.open_page('metainf', task_id=self.current_task.id)

	def bind_model_data(self, ui_element, attr):
		if not self.current_task:
			return None
		prefix = '.'.join(attr.split('.')[:-1])
		attr = attr.split('.')[-1]
		if prefix:
			prefix_obj = attrgetter(prefix)(self.current_task.data)
		else:
			prefix_obj = self.current_task.data
		ui_element.bind_value(prefix_obj, attr)

	def get_model_data(self, item):
		if not self.current_task:
			return None
		return attrgetter(item)(self.current_task.data)

	async def query_model_data(self, items: list | str):
		if not self.current_task:
			return None
		if type(items) == str:
			items = [items]

		queries = [self.current_task.data.query(item) for item in items]
		result = await asyncio.gather(*queries)
		if len(result) == 1:
			return result[0]
		return result

	def reload_model(self) -> None:
		self.current_task.data.reload()

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
