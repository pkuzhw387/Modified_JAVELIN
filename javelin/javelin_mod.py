from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model, Pmap_Model


cy = get_data(["/home/zhanghaowen/r.txt", "/home/zhanghaowen/g.txt"])
cymod = Pmap_Model(GPmodel="DRW", zydata=cy)
cymod.do_mcmc()
cymod.show_hist()

