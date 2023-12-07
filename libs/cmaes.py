import sys
sys.path.append("/usr/local/lib")
import lcmaes

stopping_criteria_dict = {'CONDITIONCOV':-15, 'TOLHISTFUN':1, 'TOLX':2, 'NOEFFECTAXIS':3, 'NOEFFECTCOOR':4, \
                          'EQUALFUNVALS':5, 'STAGNATION':6, 'AUTOMAXITER':7, 'MAXFEVALS':8, \
                          'MAXITER':9, 'FTARGET':10}
                          
class CMAES:
    def __init__(self, lambda_, seed, sigma, stopping_criteria, max_iter, outfile):
        self.lambda_ = lambda_
        self.seed = seed
        self.sigma = sigma
        self.stopping_criteria = stopping_criteria
        self.max_iter = max_iter
        self.dim = 3
        self.outfile = outfile

    def optimize(self, fitfunc, init_param):
        p  = lcmaes.make_simple_parameters(init_param, self.sigma, self.lambda_, self.seed)
        self.set_stopping_criteria(p)
        p.set_max_iter(self.max_iter)
        objfunc = lcmaes.fitfunc_pbf.from_callable(fitfunc)
        cmasols = lcmaes.pcmaes(objfunc,p)
        bcand = cmasols.best_candidate()
        bx = lcmaes.get_candidate_x(bcand)
        print("best x=",bx)
        print("distribution mean=",lcmaes.get_solution_xmean(cmasols))
        cov = lcmaes.get_solution_cov(cmasols) # numpy array
        print("cov=",cov)
        print("elapsed time=",cmasols.elapsed_time(),"ms")
        return bx
        
    def set_stopping_criteria(self, p):
        scd = stopping_criteria_dict.copy()
        for stopping_criteria in self.stopping_criteria:
            scd.pop(stopping_criteria)
            
        for s in scd:
            p.set_stopping_criteria(scd[s], False)
            

class CMAES_PWQB_LS(CMAES):
    def __init__(self, lambda_, seed, sigma, stopping_criteria, max_iter, outfile, lbounds, ubounds):
        super().__init__(lambda_, seed, sigma, stopping_criteria, max_iter, outfile)
        self.lbounds = lbounds
        self.ubounds = ubounds
        
    def optimize(self, fitfunc, init_param):
        gp = lcmaes.make_genopheno_pwqb_ls(self.lbounds, self.ubounds, self.dim)
        p  = lcmaes.make_parameters_pwqb_ls(init_param, self.sigma, gp, self.lambda_, self.seed)
        self.set_stopping_criteria(p)
        p.set_max_iter(self.max_iter)
        objfunc = lcmaes.fitfunc_pbf.from_callable(fitfunc)
        cmasols = lcmaes.pcmaes_pwqb_ls(objfunc,p)
        bcand = cmasols.best_candidate()
        bcand_pheno = lcmaes.get_best_candidate_pheno_ls(cmasols,gp)
        bx = lcmaes.get_candidate_x(bcand)
        print("best x=",bx)
        print("best x in phenotype space=",bcand_pheno)
        print("distribution mean=",lcmaes.get_solution_xmean(cmasols))
        cov = lcmaes.get_solution_cov(cmasols) # numpy array
        print("cov=",cov)
        print("elapsed time=",cmasols.elapsed_time(),"ms")
        return bcand_pheno

