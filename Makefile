PATHTOOLS=$$HOME/NiftyPET_tools
HMUDIR=$(PATHTOOLS)/mmr_hardwareumaps
DATA=Ab_PET_mMR_test

old.%: PY=$$HOME/miniconda3/envs/niftypet_old/bin/python
old.%: NINST=6d3e7b9
old.%: NIMPA=18d0980
old.%: NIPET=0376266

new.%: PY=$$HOME/miniconda3/envs/niftypet/bin/python
new.%: NINST=devel
new.%: NIMPA=devel
new.%: NIPET=devel

%.setup:
	PATHTOOLS="$(PATHTOOLS)" HMUDIR="$(HMUDIR)" $(PY) -m pip install argopt tqdm
	PATHTOOLS="$(PATHTOOLS)" HMUDIR="$(HMUDIR)" $(PY) -m pip install "git+https://github.com/NiftyPET/NInst.git@$(NINST)#egg=ninst"
	PATHTOOLS="$(PATHTOOLS)" HMUDIR="$(HMUDIR)" $(PY) -m pip install "git+https://github.com/NiftyPET/NIMPA.git@$(NIMPA)#egg=nimpa"
	PATHTOOLS="$(PATHTOOLS)" HMUDIR="$(HMUDIR)" $(PY) -m pip install "git+https://github.com/NiftyPET/NIPET.git@$(NIPET)#egg=nipet"
old.cprof:
	OMP_NUM_THREADS=1 $(PY) -m cProfile -o $@ benchmark.py $(DATA)
new.cprof:
	OMP_NUM_THREADS=1 $(PY) -m cProfile -o $@ benchmark.py --cuvec $(DATA)
