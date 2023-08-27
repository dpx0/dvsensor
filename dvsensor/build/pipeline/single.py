import asyncio


async def analyze_single_sequence(task_handler, task_options):
	# data = [{'position': new_id,
	# 	 'triplet': new_id,
	# 	 'region': new_id,
	# 	 'range': new_id,
	# 	 'percent_gc': new_id,
	# 	 'n_stop_codons': new_id,
	# 	 'off_targets': new_id}
	# 		for new_id in range(50)]

	ui_connection = task_handler.ui_connection

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
	progress_bar = ui_connection.get_element('progress_bar')
	while True:
		await asyncio.sleep(2)
		if not data:
			print("no data")
			break
		ui_connection.call('add_rows', [data.pop()])

		ui_connection.call('update_progress', step=progress_step)
		print(progress_bar.value)
		if round(progress_bar.value * 100) >= 100:
			print("FINISHED")
			ui_connection.call('set_status_finished')
			task_handler.terminate()
			break
