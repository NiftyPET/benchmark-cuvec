PATHTOOLS=$$HOME/NiftyPET_tools
HMUDIR=$(PATHTOOLS)/mmr_hardwareumaps
DATA=Ab_PET_mMR_test

old.%: PY=$$HOME/miniconda3/envs/old/bin/python
old.%: NINST=6d3e7b9
old.%: NIMPA=18d0980
old.%: NIPET=0376266

new.%: PY=$$HOME/miniconda3/envs/new/bin/python
new.%: NINST=devel
new.%: NIMPA=devel
new.%: NIPET=devel

ENV=PATHTOOLS="$(PATHTOOLS)" HMUDIR="$(HMUDIR)" OMP_NUM_THREADS=1
%.setup:
	$(ENV) $(PY) -m pip install argopt tqdm
	$(ENV) $(PY) -m pip install "git+https://github.com/NiftyPET/NInst.git@$(NINST)#egg=ninst"
	$(ENV) $(PY) -m pip install "git+https://github.com/NiftyPET/NIMPA.git@$(NIMPA)#egg=nimpa"
	$(ENV) $(PY) -m pip install "git+https://github.com/NiftyPET/NIPET.git@$(NIPET)#egg=nipet"
old.cprof:
	$(ENV) $(PY) -m cProfile -o $@ benchmark.py $(DATA)
new.cprof:
	$(ENV) $(PY) -m cProfile -o $@ benchmark.py --cuvec $(DATA)
