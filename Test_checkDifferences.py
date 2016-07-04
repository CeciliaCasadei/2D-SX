import os
gitFiles = os.listdir('./2d_mx_processing-1e4175dc69fce8c3a4d62735890660d663470954-1e4175dc69fce8c3a4d62735890660d663470954')
for gitFile in gitFiles:
	os.system('diff -u ./2d_mx_processing-1e4175dc69fce8c3a4d62735890660d663470954-1e4175dc69fce8c3a4d62735890660d663470954/%s ../2d_mx_processing/%s > ./diff_%s'%(gitFile, gitFile, gitFile))
