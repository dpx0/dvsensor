from nicegui import background_tasks
from model.taskdata import TaskData
from pipeline.single import analyze_single_sequence
import asyncio
import uuid


class TaskHandler:

	def __init__(self, controller, task_data: TaskData) -> None:
		self.controller = controller
		self.data = task_data
		self.id = str(uuid.uuid4())
		self.ui_connection = None
		self.task = None
		self.thread = None
		self.ui_connection_established = None

	def start(self, task_options) -> None:
		self.ui_connection_established = asyncio.Event()
		self.task = background_tasks.create(self._task(task_options))

	def cancel(self) -> None:
		if self.task:
			self.task.cancel()
			self.terminate()

	def terminate(self) -> None:
		self.task = None
		self.thread = None
		self.ui_connection = None
		self.ui_connection_established = None
		self.data = None

	async def _task(self, task_options) -> None:
		await self.ui_connection_established.wait()
		await self.ui_connection.page_render_complete.wait()

		if task_options['type'] == 'single':
			self.thread = asyncio.to_thread(analyze_single_sequence, self, task_options)
			await self.thread

	def connect_ui(self, ui_connection) -> None:
		if self.task:
			self.ui_connection_established.set()
			self.ui_connection = ui_connection
			self.ui_connection.add_function('cancel_task', self.ui_cancel_task)

	def ui_cancel_task(self, *args, **kwargs) -> None:
		self.cancel()
