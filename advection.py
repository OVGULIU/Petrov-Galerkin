#! /usr/bin/python3

from nutils import mesh, cli, log, function, plot, numeric, solver, _
import numpy, unittest


class MakePlots:

  def __init__(self, domain):
    self.domain = domain
    self.index = 0

  def __call__(self, ns):
    self.index += 1
    xp, up = self.domain.elem_eval([ns.x, ns.u], ischeme='bezier7', separate=True)
    with plot.PyPlot('solution', index=self.index) as plt:
      plt.mesh(xp, up)
      plt.clim(0, 1)
      plt.colorbar()


def main(
    nelems: 'number of elements' = 20,
    degree: 'polynomial degree' = 1,
    timescale: 'time scale (timestep=timescale/nelems)' = .5,
    tol: 'solver tolerance' = 1e-5,
    ndims: 'spatial dimension' = 1,
    endtime: 'end time, 0 for no end time' = 0,
    figures: 'create figures' = True,
 ):

  # construct mesh, basis
  ns = function.Namespace()
  domain, ns.x = mesh.rectilinear([numpy.linspace(0,1,nelems+1)]*ndims, periodic=range(ndims))
  ns.basis = domain.basis('discont', degree=degree)
  ns.u = 'basis_n ?lhs_n'

  # construct initial condition (centered gaussian)
  lhs0 = domain.project('exp(-?y_i ?y_i)(y_i = 5 (x_i - 0.5_i))' @ ns, onto=ns.basis, geometry=ns.x, degree=5)

  # prepare residual
  ns.f = '.5 u'
  ns.C = 1
  res = domain.integral('-basis_n,0 f' @ ns, geometry=ns.x, degree=5)
  res += domain.interfaces.integral('-[basis_n] n_0 ({f} - .5 C [u] n_0)' @ ns, geometry=ns.x, degree=5)
  inertia = domain.integral('basis_n u' @ ns, geometry=ns.x, degree=5)

  # prepare plotting
  makeplots = MakePlots(domain) if figures else lambda ns: None

  # start time stepping
  timestep = timescale/nelems
  for itime, lhs in log.enumerate('timestep', solver.impliciteuler('lhs', res, inertia, timestep, lhs0, newtontol=tol)):
    makeplots(ns(lhs=lhs))
    if endtime and itime * timestep >= endtime:
      break

  return res.eval(arguments=dict(lhs=lhs)), lhs

if __name__ == '__main__':
  cli.run(main)
