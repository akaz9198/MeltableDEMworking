#!BPY

__doc__ = """
AnimateGrains.py

This script animates some list of grains as follows:
1) Delete (unlink) and replacement of default objects in a scene
2) Addition of spheres from a data file from Pierre's DEM code
2) Using python to render a movie or still of the spheres in Blender

This script is executed at the command line by:
>blender -P AnimateGrains.py -- [options]
"""
__author__ = "benjy and tom"
__version__ = "0.7 03/05/2010"

# Add a licence here if you wish to re-distribute, we recommend the GPL

from Blender import Ipo, Mathutils
import bpy, BPyMessages, math
from math import *
import sys
import optparse

def makeParticleIpo(object, firstFrame, numberOfFrames, frameStep, grains, pt, chain):
	# Create an new Ipo Curve of name myIpo and type Object
	myIpo = bpy.data.ipos.new(str(pt), 'Object')
	
	# Create LocX, LocY, and LocZ Ipo curves in our new Curve Object
	# and store them so we can access them later
	myIpo_x = myIpo.addCurve('LocX')
	myIpo_y = myIpo.addCurve('LocY')
	myIpo_z = myIpo.addCurve('LocZ')
	if chain != 'off':
		myIpo_layer = myIpo.addCurve('Layer')
	
	# This Calculates the End Frame for use in an xrange() expression
	endFrame = firstFrame + (numberOfFrames * frameStep) + frameStep
	
	for frame in xrange(firstFrame, endFrame, frameStep):
		# check if there is a force chain
		if chain != 'off':
			chain_length = 2 # minimum chain length to plot
			if chain[frame-1,pt] > chain_length:
				ipoValue_layer = 1
			else:
				ipoValue_layer = 2

		ipoValue_x = grains[pt,frame-1,0]
		ipoValue_y = grains[pt,frame-1,1]
		ipoValue_z = grains[pt,frame-1,2]
		
		# Append to the Ipo curve at location frame, with the value ipoValue_x
		# Note that we should pass the append function a tuple or a BezTriple
		myIpo_x.append((frame, ipoValue_x))
		myIpo_y.append((frame, ipoValue_y))
		myIpo_z.append((frame, ipoValue_z))

		if chain != 'off':
			myIpo_layer.append((frame, ipoValue_layer))	
			
	# Link our new Ipo Curve to the passed object
	object.setIpo(myIpo)
	
def makeCameraObjectIpo(object, firstFrame, numberOfFrames, frameStep):
	# Create an new Ipo Curve of name myIpo and type Object
	myIpo = bpy.data.ipos.new('myIpo', 'Object')
	
	# Create LocX, LocY, and LocZ Ipo curves in our new Curve Object
	# and store them so we can access them later
	myIpo_x = myIpo.addCurve('LocX')
	myIpo_y = myIpo.addCurve('LocY')
	myIpo_z = myIpo.addCurve('LocZ')
	
	# This Calculates the End Frame for use in an xrange() expression
	endFrame = firstFrame + (numberOfFrames * frameStep) + frameStep
		
	for frame in xrange(firstFrame, endFrame, frameStep):
		t=float(frame)*math.pi/500 # one rotation every 500 frames
		ipoValue_x = 20*math.cos(2*t)
		ipoValue_y = t-float(endFrame)/20 # spiral in y direction
		ipoValue_z = 20*math.sin(2*t)
		
		# Append to the Ipo curve at location frame, with the value ipoValue_x
		# Note that we should pass the append function a tuple or a BezTriple
		myIpo_x.append((frame, ipoValue_x))
		myIpo_y.append((frame, ipoValue_y))
		myIpo_z.append((frame, ipoValue_z))	


	# Link our new Ipo Curve to the passed object
	object.setIpo(myIpo)
		
def makeCameraLensIpo(object, firstFrame, numberOfFrames, frameStep):
	# Create an new Ipo Curve of name myIpo and type Object
	myIpo = bpy.data.ipos.new('myIpo', 'Camera')
	
	# Create LocX, LocY, and LocZ Ipo curves in our new Curve Object
	# and store them so we can access them later
	myIpo_Lens = myIpo.addCurve('Lens')
	
	# This Calculates the End Frame for use in an xrange() expression
	endFrame = firstFrame + (numberOfFrames * frameStep) + frameStep
		
	for frame in xrange(firstFrame, endFrame, frameStep):
		LensMax = 40. # how zoomed out at start
		LensMin = 30. # how zoomed in at end
		ipoValue_Lens = float(frame)*(LensMin - LensMax)/(endFrame - 1) + LensMax
		myIpo_Lens.append((frame, ipoValue_Lens))
	
	# Link our new Ipo Curve to the passed object
	object.setIpo(myIpo)
	
def buildParticles(grains,particles,q):
	import csv
	from Blender import Camera, Lamp, Material, Mesh, Object, Constraint
		
	scene = bpy.data.scenes.active
    ##############################################################
    # Get rid of all object from default scene
	print 'Removing objects in current scene'
	for ob in scene.objects:
		if ((cmp(ob.getName(),'Lamp')==0), (cmp(ob.getName(),'Cube')==0), (cmp(ob.getName(),'Camera')==0)):
			scene.objects.unlink(ob)
			print str(ob) + 'removed'
	print 'Adding new objects'
    ##############################################################
    # add mesh at origin to focus camera on
    #
	me = Mesh.Primitives.UVsphere(3,3,0)
	origin = scene.objects.new(me,'Origin')
    ##############################################################
    # add a camera and set it up
    #
	camdata = Camera.New('ortho')   # create new ortho camera data
#	camdata.scale = 25.0             # set scale value for ortho view
	cam = scene.objects.new(camdata)   # add a new camera object from the data
	scene.objects.camera = cam    # make this camera the active camera
    
    # This makes the camera follow the origin
	ob = Object.Get('Camera')
	const = ob.constraints.append(Constraint.Type.TRACKTO)
	const[Constraint.Settings.TARGET] = origin
	const.__setitem__(Constraint.Settings.TRACK, Constraint.Settings.TRACKNEGZ)
	const.__setitem__(Constraint.Settings.UP, Constraint.Settings.UPY)

    ##############################################################
    # add a lamp and set it up
    #
	lampdata = Lamp.New()
	lampdata.setEnergy(1.0)
#	lampdata.setType('Sun')
	lamp = scene.objects.new(lampdata)
	lamp.setLocation(10.0,-10.0,10.0)
	lamp.setEuler(0.0,0.0,0.0)
	##############################################################
    # get particle data
    #
	for pt in range(particles):
		me = Mesh.Primitives.UVsphere(q,q,grains[pt,0,3])   # create a new sphere of (segments,rings,diameter)
		ob = scene.objects.new(me,'Grain' + str(pt))        # add a new mesh-type object to the scene
		ob.setLocation (grains[pt,0,0], grains[pt,0,1], grains[pt,0,2])       # position the object in the scene
		# smooth the mesh
		faces = me.faces            # Get the faces.
		for f in faces:             # Loop through them.
			f.smooth = 1            # Smooth each face.
		# set the material
		mat = Material.New('Mat') # create a new Material called 'Mat'
		mat.rgbCol = [grains[pt,0,4], grains[pt,0,5], grains[pt,0,6]] # change the color of the sphere
		mat.setAlpha(0.8)                     # mat.alpha = 0.2 -- almost transparent
		mat.emit = 0.5                        # equivalent to mat.setEmit(0.8)
		mat.mode |= Material.Modes.ZTRANSP    # turn on Z-Buffer transparency
		mat.setAdd(1) 
		me.materials += [mat]

def getParticleData(particles, file_end, file_start, data_file_dir):
	##############################################################
	# loop through data files
	#
	from numpy import zeros
	import csv
	grains = zeros((particles,file_end-file_start+1,8))
	print 'Loading position files'
	for data in range(file_start,file_end+1):
		pos	= csv.reader(open(data_file_dir + 'pos_'+ str(data)), delimiter='\t') # open data
		pt = 0
		for line in pos:
			grains[pt,data-file_start,0] = float(line[0]) # find x value
			grains[pt,data-file_start,1] = float(line[1]) # find y value
			grains[pt,data-file_start,2] = float(line[2]) # find z value
			grains[pt,data-file_start,3] = 2*float(line[9]) # find diameter
			if float(line[11]) == 0:
				grains[pt,data-file_start,4] = 1 # find r color intensity
				grains[pt,data-file_start,5] = 1 # find g color intensity
			else:
				grains[pt,data-file_start,6] = 1 # find b color intensity
			grains[pt,data-file_start,7] = float(line[7]) # find stress???
			pt += 1
	return grains
	
def getContactData(data,particles,file_end,file_start,data_file_dir):
	import csv
	cont = csv.reader(open(data_file_dir + 'cont_'+ str(data)), delimiter='\t') # open data
	A=[]
	B=[]
	for line in cont:
		A.append(line[0])
		B.append(line[1])
	return A, B

# call this for each timestep to work out contacts for particle pt
def neighbours(pt, A, B, threshold, grains, time):
	neighbour_list = []
	for a, b in zip(A, B):
		if int(a) == int(pt) and grains[int(b),time,7] >= threshold:
			neighbour_list.append(int(b))
		elif int(b) == int(pt) and grains[int(a),time,7] >= threshold:
			neighbour_list.append(int(a))
	return neighbour_list
		
# returns length of force chain including particle pt
def force_chain(A, B, threshold, grains, time, particles, chain):
	for pt in range(particles):
		neighbours_list = neighbours(pt, A, B, threshold, grains, time)
		if neighbours_list:
			for i in neighbours_list:
				# find if the neighbours have neighbours. Do this for i<=chain_length
				#if neighbours(int(neighbours_list[i]), threshold):
				chain[time,pt] += 1
	return chain

def main(argv):

	# get the args passed to blender after "--", all of which are ignored by blender specifically
	# so python may receive its own arguments

	if '--' not in argv:
		argv = [] # as if no args are passed
	else:
		argv = argv[argv.index('--')+1: ] # get all args after "--"
	
	# When --help or no args are given, print this help
	usage_text =  'Run blender with this script:'
	usage_text += '  blender -P  AnimateGrains.py  -- [options]'

	parser = optparse.OptionParser(usage = usage_text)

	# Possible types are: string, int, long, choice, float and complex.
	parser.add_option('-d', '--directory', dest='data_file_dir', help='This is the directory of the data files', metavar='string')
	parser.add_option('-s', '--start', dest='file_start', help='This is the first file of data positions', metavar='int')
	parser.add_option('-f', '--finish', dest='file_end', help='This is the last file of data positions', metavar='int')
	parser.add_option('-e', '--export', dest='export_dir', help='This is the directory to export renders to', metavar='string')
	parser.add_option('-q', '--quality', dest='q', help='This is the quality of the  individual spheres', metavar='int')
	parser.add_option('-t', '--particles', dest='particles', help='This is the file containing the contacts', metavar='int')
	parser.add_option('-r', '--render', dest='render', help='Render an image to the specified path', metavar='string')
	parser.add_option('-c', '--chain', dest='chain_mode', help='Only show particles in a force chain', metavar='string')

	options, args = parser.parse_args(argv) # we don't use the args

	##############################################################
	# force a quit if required input not given
	#
	if not argv:
		parser.print_help()
		sys.exit()
	if not options.file_start:
		print 'Error: --start="some file" argument not given, aborting.'
		parser.print_help()
		sys.exit()
	if not options.file_end:
		print 'Error: --finish="some file" argument not given, aborting.'
		parser.print_help()
		sys.exit()
	if not options.particles:
		print 'Error: --particles="some number" argument not given, aborting.'
		parser.print_help()
		sys.exit()
	if not options.data_file_dir:
		print 'Error: --data_file_dir="some file" argument not given, aborting.'
		parser.print_help()
		sys.exit()
	##############################################################
	# inputs
	#
	data_file_dir = options.data_file_dir
	file_start = int(options.file_start)
	file_end = int(options.file_end)
	particles = int(options.particles)
	
	if options.q:
		q = int(options.q)
	else:
		q = 10;
	if options.render:
		render = options.render
	else:
		render = 'no'
	if options.export_dir:
		export_dir = options.export_dir
	else:
		export_dir = '//renders/'
	if options.chain_mode:
		chain_mode = options.chain_mode
	else:
		chain_mode = 'off' # off or on

	threshold = 0.001
	
	from Blender.Scene import Render
	# Get the active scene, since there can be multiple ones
	scene = bpy.data.scenes.active
	context = scene.getRenderingContext()
	context.extensions = True
	context.sizePreset(Render.PREVIEW)
	context.imageType = Render.AVIRAW
	
	# Get data from source files
	grains = getParticleData(particles,file_end,file_start,data_file_dir)
	if chain_mode != 'off':
		print 'Loading contact files'
		from numpy import zeros
		chain = zeros((file_end-file_start+1,particles))
		for data in range(file_start,file_end+1):
			A, B = getContactData(data,particles,file_end,file_start,data_file_dir)
			chain = force_chain(A, B, threshold, grains, data-file_start, particles,chain)
	else:
		chain = 'off'
				
	# Build initial particles and materials
	buildParticles(grains,particles,q)
	
	# Create IpoCurve for each particle
	from Blender import Object, Camera
	print 'Creating IpoCurves for particles'
	context.sFrame = 1
	context.eFrame = file_end - file_start
	for pt in range(particles):
		object = Object.Get('Grain' + str(pt))
		makeParticleIpo(object, scene.render.sFrame, scene.render.eFrame, 1, grains, pt, chain)
	
	# Create IpoCurve for camera object
	object = Object.Get('Camera')
	print 'Creating IpoCurve for camera object'
	makeCameraObjectIpo(object, scene.render.sFrame, scene.render.eFrame, 1)

	# Create IpoCurve for camera camera!
	object = Camera.Get('Camera')
	print 'Creating IpoCurve for camera lens'
	makeCameraLensIpo(object, scene.render.sFrame, scene.render.eFrame, 1)

	
	# Render movie
	if render == 'yes' or render == 'y':
		print 'Rendering movie'
		context.renderPath = export_dir
		context.renderAnim()
	else:
		print 'Not rendering a movie'
		

if __name__ == '__main__':
	main(sys.argv)
