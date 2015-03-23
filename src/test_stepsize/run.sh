# Meant for a machine with 3 GPUs, e.g. jak
export CUDA_VISIBLE_DEVICES=0
python production_amber.py 0.5&

sleep 2
export CUDA_VISIBLE_DEVICES=1
python production_amber.py 1.0&

sleep 2
export CUDA_VISIBLE_DEVICES=2
python production_amber.py 2.0&
python production_amber.py 1.5&
