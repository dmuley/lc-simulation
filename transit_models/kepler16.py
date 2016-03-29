if __name__ == '__main__':
    if __package__ is None:
        import sys
        from os import path
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) );
	print sys.path[-1];
        from gravsims import basicKepler as bk;
	from integrators import polarIntegrator3 as pInt;
    else:
	from ..gravsims import basicKepler as bk;
	from ..integrators import polarIntegrator3 as pInt;


