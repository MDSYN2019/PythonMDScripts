import sys
import os
import click
import subprocess
from subprocess import call

# Example: -w output 1lyl.pdb AC 12as.pdb AB

#

IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"


# Calling perl for ialign

#subprocess.call(IalignPath)

def callialign(folder, pdb1, chain1, pdb2, chain2, Ialignpath = None): # Default None for Ialignpath for now
	"""
	
	"""
	IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"
	try: 
		p = subprocess.Popen([IalignPath, "-w", folder, pdb1, chain1, pdb2, chain2], stdout = subprocess.PIPE) # Ok this works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot read ialign properly. Please check the input file/chain/pdb files, or check the path to the ialign binary")
	else:
		print ("Running ialign sucessfully..")

#@click.command()
#@click.option('--ialignbin', default = IalignPath, help = "The path to where ialign.pl is located on your computer")
#@click.argument('ialignpath')
#def call_ialign(ialignpath):
#	""" Simple program that calls ialign """
#	click.echo("Calling from path {}".format(str(ialignbin)))
#	subprocess.call(ialignbin)


