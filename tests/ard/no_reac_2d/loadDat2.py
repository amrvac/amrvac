import yt 
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
plt.rc('font', family='serif', size=20)
plt.rcParams['mathtext.fontset'] = "stix"
plt.ion()
# Load the dataset.



units={}
units["m1"]=r"$kg/m^2/s$"
units["rho1"]=r"$kg/m^3$"
units["rho"]=r"$kg/m^3$"
units["v1"]=r"$m/s$"
units["bx1"]=r"$T$"

labels={}
labels["m1"]=r"$\rho v_z$"
labels["rho1"]=r"$\rho_1$"
labels["rho"]=r"$\rho$"
labels["v1"]=r"$v_z$"
labels["bx1"]=r"$B_{\rm x1}$"

from os.path import join


dataPert = None


def v_cn(field,data):
  return np.sqrt(((data["m_c1"]/data["rho_c"] - data["m_n1"]/data["rho_n"]).to_ndarray())**2+
         ((data["m_c2"]/data["rho_c"] - data["m_n2"]/data["rho_n"]).to_ndarray())**2)

def vx_cn(field,data):
  return data["m_c2"]/data["rho_c"] - data["m_n2"]/data["rho_n"]
def vz_cn(field,data):
  return data["m_c2"]/data["rho_c"] - data["m_n2"]/data["rho_n"]


def v_n(field,data):
  return np.sqrt(((data["m_n1"]/data["rho_n"]).to_ndarray())**2+
         ((data["m_n2"]/data["rho_n"]).to_ndarray())**2)

def v_n1(field,data):
  return  (data["m_n2"]/data["rho_n"]).to_ndarray()

def v_n2(field,data):
  return  (data["m_n2"]/data["rho_n"]).to_ndarray()


def v_c(field,data):
  return np.sqrt(((data["m_c1"]/data["rho_c"]).to_ndarray())**2+
         ((data["m_c2"]/data["rho_c"]).to_ndarray())**2)


def v_c1(field,data):
  return  (data["m_c1"]/data["rho_c"]).to_ndarray() / (
   #np.sqrt((data["b1"].to_ndarray()**2 + data["b2"].to_ndarray()**2 + data["b3"].to_ndarray()**2 )/data["rho_c"]) )
   np.sqrt((data["b1"].to_ndarray()**2 + data["b2"].to_ndarray()**2  )/data["rho_c"]) )

def v_c2(field,data):
  return  (data["m_c2"]/data["rho_c"]).to_ndarray()

def rho_c1(field,data):
  return  (data["rho_c"]).to_ndarray()

def e_c1(field,data):
  return  (data["p_c"]).to_ndarray()


units={}
units["vz_cn"] = 'code_length/code_time'
units["vz_cn"] = ''
units["v_cn"] = ''

from matplotlib._cm import cubehelix as _cubehelix
name=""
#s=-0.3 #blue
#s=-1.07 #green
#s=0.9 #red
#s=0.0 #violet/blue

cm1={}
cm1["u"] = "binary"
cm1["v"] = "binary"
cm1["w"] = "binary"

dir1="output"

dataset="u"

listFiles = range(0,10000)
xprobmin1     = 0e0
xprobmax1     = 64e0
xprobmin2     = 0e0
xprobmax2     = 64e0

center = (0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),0)
width = (xprobmax1-xprobmin1,xprobmax2-xprobmin2,0)


from os.path import exists
for k in listFiles:
  fn =  "%s/pure_adv_%04d.dat" % (dir1,k)
  if(not exists(fn)):
    continue
  ds = yt.load(fn)
  print("FIELDS saved ", ds.field_list)
  data = ds.all_data()
  print("FIELDS ",ds.field_list)
  if(not ('amrvac',dataset) in ds.field_list):
    #ds.add_field(('gas', dataset), function=locals()[dataset], units=units[dataset], sampling_type='cell')
    ds.add_field(('gas', dataset), function=locals()[dataset], units='', sampling_type='cell')


  print

  #clipPlane = ds.slice("y", limSup)

  
  print(k, "***HASNAN***", np.isnan(np.sum(data[dataset])))
  print("---MINMAX---", np.min(data[dataset]), np.max(data[dataset])  )

  p = yt.ProjectionPlot(ds, "z",dataset, origin="native", center=center, width=width,method="sum")


  #p.set_axes_unit("m")
  p.hide_colorbar()

  p.set_cmap(dataset,cm1[dataset])
  #p.set_zlim(dataset, 1.0)
  p.set_log(dataset, False)

  print("DIR P ", dir(p))
  #p.annotate_grids()
  #p.annotate_streamlines('b1','b2' ,plot_args={"linewidth": 0.25,"arrowsize":0.25, "color":"blue"},density=0.75)
  #p.annotate_quiver(('amrvac',"v_c1"), ('amrvac',"v_c2"), factor=16,
  #                plot_args={"color": "orange"})
  p.show_colorbar()

  p._setup_plots()
  cb = p.plots[dataset].cb
  cb.set_label(p._colorbar_label[dataset],rotation=-90,labelpad=20)
  ax = p.plots[dataset].axes  
  ax.set_xlabel("z")
  ax.set_ylabel("x", rotation=90)


  ax.invert_yaxis() 
  #lx = 0.08
  lx = 0.92
  ly = 0.92
  #ly = 0.06
  color="black"
  #color="lime"
  #color="lime"
  #color="tab:pink"
  #unit_time = 155.64911177469682
  ax.text(lx, ly,  "t=%.1f"%(ds["time"]),size=14, color=color,ha="center", va="center",transform=ax.transAxes, rotation=-90 )
  p.save(join(dir1,("%s-%04d.png"% (dataset,k))))
  #p.save(join("~/IMG",("%s-%04d.png"% (dataset,k))))
