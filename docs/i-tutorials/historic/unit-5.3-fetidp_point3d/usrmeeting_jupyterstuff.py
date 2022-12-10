uid_prefix = 'STDID'
import time, sys, subprocess
import ipyparallel as ipp

def start_cluster(n, user_id='', profile='ngsolve'):
    user_id = uid_prefix+user_id
    if profile!='ngsolve' and uid_prefix not in ['kogler', 'lukas']:
        print('Please use profile ngsolve. ngsolve2 spams the login node with MPI jobs.')
        quit()
    ret = subprocess.call(['ipcluster','start','-n='+str(n),'--profile='+profile,'--cluster-id='+str(user_id),'--daemonize=True'])
    if ret != 0:
        print('could not start cluster, (probably already/still running)')

def stop_cluster(user_id='', profile='ngsolve'):
    user_id = uid_prefix+user_id
    if profile!='ngsolve' and uid_prefix not in ['kogler', 'lukas']:
        print('Please use profile ngsolve. ngsolve2 spams the login node with MPI jobs.')
        quit()
    ret = subprocess.call(['ipcluster','stop','--profile='+profile,'--cluster-id='+str(user_id),'--daemonize=True'])
    time.sleep(3) 
    # if ret != 0:
    #     print('could not stop cluster, (probably not running)')
    # else:
    #     print('successfully stopped cluster!')

def connect_cluster(user_id='', profile='ngsolve'):
    user_id = uid_prefix+user_id
    if profile!='ngsolve' and uid_prefix not in ['kogler', 'lukas']:
        print('Please use profile ngsolve. ngsolve2 spams the login node with MPI jobs.')
        quit()
    ippclient = ipp.Client(profile=profile, cluster_id=user_id, block=True)
    ippclient.activate()
    for k in range(300):
        try:
            ippclient[:]
        except:
            sys.stdout.write('\rconnecting ... try:'+str(k))
            sys.stdout.flush()
            time.sleep(1)
            continue
        else:
            sys.stdout.write('\rconnecting ... try:'+str(k)+str(' succeeded!'))
            sys.stdout.flush()
            time.sleep(3)
            return
    raise Exception('Could not connect, try again!')
