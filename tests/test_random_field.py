import numpy as np
import ucla_geotech_tools.random_field as rf

def test_random_field1():
    z = np.random.rand(1000)
    x = np.random.rand(1000)
    theta_z = 0.5
    mu = 0
    sigma = 1
    model = 0
    field1a = rf.get_random_field(x=x,z=z,theta_z=theta_z,mu=mu,sigma=sigma,model=model)
    model = 1
    field1b = rf.get_random_field(x=x,z=z,theta_z=theta_z,mu=mu,sigma=sigma,model=model)
    
def test_random_field2():
    z = np.random.rand(1000)
    x = np.random.rand(1000)*10
    theta_z = 0.5
    theta_x = 10*theta_z
    mu = 0
    sigma = 1
    model = 0
    field2a = rf.get_random_field(x=x,z=z,theta_z=theta_z,mu=mu,sigma=sigma,model=model,theta_x=theta_x)
    field2b = rf.get_random_field(x=x,z=z,theta_z=theta_z,mu=mu,sigma=sigma,model=model,theta_x=theta_x,Dhchol=0.1)