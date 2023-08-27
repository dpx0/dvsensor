from nicegui import background_tasks
from model.rnadata import RNAData
from pipeline.single import analyze_single_sequence
import asyncio
import uuid


class TaskHandler:

	def __init__(self, controller, rna_data: RNAData) -> None:
		self.controller = controller
		self.rna_data = rna_data
		self.id = str(uuid.uuid4())
		self.ui_connection = None
		self.async_task = None
		self.ui_connection_flag = None

	@property
	def model(self):
		return self.rna_data

	def start(self, task_options) -> None:
		print(f"starting task {self.id} with options {task_options}")
		self.ui_connection_flag = asyncio.Event()
		self.async_task = background_tasks.create(self._task(task_options))

	def cancel(self) -> None:
		if self.async_task:
			self.async_task.cancel()
			self.terminate()

	def terminate(self) -> None:
		self.async_task = None
		self.ui_connection = None
		self.ui_connection_flag = None

	async def _task(self, task_options) -> None:
		await self.ui_connection_flag.wait()
		print("starting task...")
		print(self.ui_connection)

		if task_options['type'] == 'single':
			await analyze_single_sequence(self, task_options)

	def reload_model(self) -> None:
		original_record = self.rna_data.original_sequence_record
		self.rna_data = RNAData(original_record)

	def connect_ui(self, ui_connection) -> None:
		if self.async_task:
			self.ui_connection_flag.set()
			self.ui_connection = ui_connection
			self.ui_connection.add_function('cancel_task', self.ui_cancel_task)

	def ui_cancel_task(self, *args, **kwargs) -> None:
		print("cancelling task...")
		self.cancel()

