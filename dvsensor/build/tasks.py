from nicegui import background_tasks
import asyncio
from model.rnadata import RNAData
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

		# data = [{'position': new_id,
		# 	 'triplet': new_id,
		# 	 'region': new_id,
		# 	 'range': new_id,
		# 	 'percent_gc': new_id,
		# 	 'n_stop_codons': new_id,
		# 	 'off_targets': new_id}
		# 		for new_id in range(50)]

		data = [
			{'position': '7127',
			 'triplet': 'CCA',
			 'region': '3UTR',
			 'range': '7079-7178',
			 'percent_gc': '52,8%',
			 'n_stop_codons': '0',
			 'off_targets': '0'},
			{'position': '5788',
			 'triplet': 'CCA',
			 'region': '3UTR',
			 'range': '5740-5839',
			 'percent_gc': '51,2%',
			 'n_stop_codons': '0',
			 'off_targets': '0'},
		]

		progress_step = 1.0 / len(data)
		progress_bar = self.ui_connection.get_element('progress_bar')
		while True:
			await asyncio.sleep(2)
			if not data:
				print("no data")
				break
			self.ui_connection.call('add_rows', [data.pop()])

			self.ui_connection.call('update_progress', step=progress_step)
			print(progress_bar.value)
			if round(progress_bar.value * 100) >= 100:
				print("FINISHED")
				self.ui_connection.call('set_status_finished')
				self.terminate()
				break

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

