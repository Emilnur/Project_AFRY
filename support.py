# === Support Functions ===
import numpy as np
import os
import re
import glob
import sys
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import openOdb
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch

# ======= DOE Generation ======= #
def create_output_structure():
    cwd = os.getcwd()
    # Filter only directories that start with 'doe_'
    doe_folders = sorted(
        [d for d in os.listdir(cwd)
         if os.path.isdir(os.path.join(cwd, d)) and d.startswith('doe_')]
    )
    next_doe_index = len(doe_folders) + 1
    doe_name = f"doe_{next_doe_index}"
    doe_path = os.path.join(cwd, doe_name)
    os.makedirs(doe_path)
    return doe_path, next_doe_index

def create_run_folder(base_path, run_index):
    run_name = f"run_{run_index:04d}"
    run_path = os.path.join(base_path, run_name)
    os.makedirs(run_path)
    return run_path

def load_doe(filepath):
    designs = []
    try:
        with open(filepath, newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                thickness = np.array([
                    float(row["Thickness_Vessel"]),
                    float(row["Thickness_Pad"]),
                    float(row["Thickness_Nozzle"])
                ])
                VR = float(row["VR"])
                NR = float(row["NR"])
                H = float(row["H"])
                LN = float(row["LN"])
                RP = float(row["RP"])
                forces = np.array([
                    float(row["Force_X"]),
                    float(row["Force_Y"]),
                    float(row["Force_Z"])
                ])
                moments = np.array([
                    float(row["Moment_X"]),
                    float(row["Moment_Y"]),
                    float(row["Moment_Z"])
                ])
                P = float(row["Pressure"])
    
                designs.append((VR, NR, H, LN, RP, thickness, forces, moments, P))
        return designs
    except:
        pass

def del_abq_temp():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # Fallback for interactive environments
        script_dir = os.getcwd()

    # Patterns to match
    patterns = [
        "abaqus.rpy",         # plain file
        "abaqus.rpy.*",       # numbered files
        "abaqus*.rec"         # rec files like abaqus1.rec, abaqus2.rec, etc.
    ]

    # Delete matching files
    for pattern in patterns:
        for file_path in glob.glob(os.path.join(script_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

    # Delete __pycache__ folder
    pycache_dir = os.path.join(script_dir, "__pycache__")
    if os.path.exists(pycache_dir) and os.path.isdir(pycache_dir):
        try:
            shutil.rmtree(pycache_dir)
            print(f"Deleted folder: {pycache_dir}")
        except Exception as e:
            print(f"Failed to delete folder {pycache_dir}: {e}")

def select_next_doe_csv():
    current_dir = os.getcwd()

    # Find all doe_*.csv files in the current directory
    csv_files = sorted(
        [f for f in os.listdir(current_dir) if re.match(r"doe_\d+\.csv$", f)],
        key=lambda x: int(re.search(r"\d+", x).group())
    )

    if not csv_files:
        print("No design of experiments file (doe_*.csv) found in the current directory.")
        return None

    # Find all folders named doe_*
    doe_folders = {
        f for f in os.listdir(current_dir)
        if os.path.isdir(f) and re.match(r"doe_\d+$", f)
    }

    # Look for the first csv file without a matching folder
    for csv_file in csv_files:
        folder_name = csv_file.rsplit('.', 1)[0]
        if folder_name not in doe_folders:
            print(f"Selected DOE file: {csv_file}")
            return csv_file

    print("No design of experiments left to run — all corresponding folders exist.")
    return None

# ======= Pre-Processing ======= #
def validate_geom(VR, NR, H, LN, RP, thickness):
    """
    Validates and sanitizes geometry input parameters.
    """
    args = locals()
    corrected_args = {}

    for name, value in args.items():
        if isinstance(value, (int, float)):
            corrected_args[name] = abs(value)
        elif isinstance(value, np.ndarray):
            corrected_args[name] = np.abs(value)
        else:
            raise TypeError(f"Unsupported input type for '{name}': {type(value)}")

    VR = corrected_args['VR']
    NR = corrected_args['NR']
    H = corrected_args['H']
    LN = corrected_args['LN']
    RP = corrected_args['RP']
    thickness = corrected_args['thickness']

    t_ves, t_pad, t_noz = thickness

    xsi_N = 2.5 * np.sqrt(NR * t_noz)
    RN = NR - t_pad/2

    if xsi_N + RN > H / 2:
        raise ValueError("Nozzle radius too large in comparison with vessel height.")

    return VR, NR, H, LN, RP, thickness, xsi_N

def Parts(mymodel, VR, NR, H, LN, RP, thickness):
    """
    Creates the geometry of the pressure vessel with nozzle and pad.
    Returns the flat base area for force calculation.
    INPUT:
    VR        - Inner radius of the pressure vessel
    NR        - Outer radius of the nozzle
    LN        - Nozzle length
    thickness - row array with the thickness of the vessel, reinforcement pad and nozzle

    OUTPUT:
    flatA     - Area of the vessel's flat base, used to compute the upper edge force   
    """
    
    # Validate geometry and unpack correct values
    VR, NR, H, LN, RP, thickness, xsi_N = validate_geom(VR, NR, H, LN, RP, thickness)

    t_ves, t_pad, t_noz = thickness
    RV = VR + t_ves/2           # midsurface radius of pressure vessel 
    RN = NR - t_pad/2           # midsurface radius of nozzle             

    #### Create Pressure Vessel Geometry ###

    # Create sketch
    s = mymodel.ConstrainedSketch(name='__profile__', 
        sheetSize=5000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)

    # Define base circle
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(RV, 0.0))

    p = mymodel.Part(name='Pressure_Vessel', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mymodel.parts['Pressure_Vessel']

    # Extrude by given Height Value
    p.BaseShellExtrude(sketch=s, depth=H)
    s.unsetPrimaryObject()

    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mymodel.sketches['__profile__']

    # Retrieves edges of the open cylinder and closes it
    e1 = p.edges
    p.CoverEdges(edgeList = e1[0:2], tryAnalytical=False)


    ### Create Nozzle Geometry ###
    # Create datum plane offset by the vessel radius + nozzle length
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=RV+LN)

    f1, e, d1 = p.faces, p.edges, p.datums

    # Define position of the center of the datum plane
    t = p.MakeSketchTransform(sketchPlane=d1[3], sketchUpEdge=e[1], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, RV+LN, 
        H/2))
    s = mymodel.ConstrainedSketch(name='__profile__', 
        sheetSize=13869.31, gridSpacing=346.73, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)

    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    # Nozzle radius as input to the circle creation
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(RN, 0.0))

    f = p.faces
    f1, e1, d2 = p.faces, p.edges, p.datums

    # Extrude nozzle circle up to pressure vesel face
    p.ShellExtrude(sketchPlane=d2[3], sketchUpEdge=e1[1], upToFace=f1[2], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, 
        flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    del mymodel.sketches['__profile__']

    f = p.faces
    p.RemoveFaces(faceList = f[4:5], deleteCells=False)

    ## Project Pad ###
    f1, e, d1 = p.faces, p.edges, p.datums

    # Define position of the center of the datum plane
    t = p.MakeSketchTransform(sketchPlane=d1[3], sketchUpEdge=e[1], 
        sketchPlaneSide=SIDE1, origin=(0.0, RV+LN, H/2))
    s1 = mymodel.ConstrainedSketch(name='__profile__', 
        sheetSize=16026.2, gridSpacing=400.65, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)

    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(RP, 0.0))

    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#2 ]', ), )
    f, e1, d2 = p.faces, p.edges, p.datums
    p.PartitionFaceBySketchDistance(sketchPlane=d2[3], sketchUpEdge=e1[1], 
        faces=pickedFaces, sketchPlaneSide=SIDE1, sketch=s1, distance=2*LN)
    s1.unsetPrimaryObject()
    del mymodel.sketches['__profile__']

    ### Sectioning for Meshing Purposes
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)

    d = p.datums
    p.DatumPlaneByRotation(plane=d[7], axis=d[8], angle=45.0)
    d1 = p.datums
    p.DatumPlaneByRotation(plane=d1[7], axis=d1[8], angle=-45.0)

    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#19 ]', ), )
    d = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d[9], faces=pickedFaces)

    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#cf ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[10], faces=pickedFaces)

    v, e = p.vertices, p.edges
    p.DatumPlaneByThreePoints(point1=v[8], point3=v[1], point2=p.InterestingPoint(
        edge=e[24], rule=CENTER))

    v1, e1 = p.vertices, p.edges
    p.DatumPlaneByThreePoints(point1=v1[5], point3=v1[9], 
        point2=p.InterestingPoint(edge=e1[24], rule=CENTER))

    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#d00 ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[14], faces=pickedFaces)

    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#6807 ]', ), )
    d = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d[13], faces=pickedFaces)

    f = p.faces
    p.RepairFaceNormals(faceList = f[0:1]+f[5:6]+f[8:9]+f[19:20])
    p.RepairFaceNormals(faceList = f[0:-1])

    f = p.faces

    # Test if face is cylindrical or not - this excludes the flat top and bottom surfaces of the pressure vessel
    pickedFaces = []
    flat_area = []
    for surf in f:
        try:
            r = surf.getRadius()
        except:
            pickedFaces.append(surf)
            area = surf.getSize()
            flat_area.append(area)
            continue

    flatA = sum(flat_area)/2
    p.RemoveFaces(faceList = pickedFaces, deleteCells=False)

    ## Create copy of part and keep pad only
    p = mymodel.Part(name='Pad', 
        objectToCopy=mymodel.parts['Pressure_Vessel'])
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p2 = mymodel.parts['Pad']
    f2 = p2.faces
    p2.RemoveFaces(faceList = f[0:1]+f[2:4]+f[5:7]+f[8:14], deleteCells=False)
    
    # Fix pad edges
    e1 = p.edges
    p.MergeEdges(edgeList = e1[10:12], extendSelection=False)
    p = mymodel.parts['Pad']
    e = p.edges
    p.MergeEdges(edgeList = e[11:13], extendSelection=False)
    mymodel.parts['Pad'].checkGeometry()

    p = mymodel.parts['Pad']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)

    # Create Sets for Pad Part
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#8a8 ]', ), )
    p.Set(edges=edges, name='Edge_Pad_Outer')

    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#642 ]', ), )
    p.Set(edges=edges, name='Edge_Pad_Inner')

    f = p.faces
    p.Set(faces=f, name='Surf_Pad_Section')

    # Fix vessel faces and edges
    p = mymodel.parts['Pressure_Vessel']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)

    e = p.edges
    p.RemoveRedundantEntities(edgeList = e[7:8]+e[13:14]+e[21:22]+e[25:29]+\
        e[32:34])

    e = p.edges
    p.MergeEdges(edgeList = e[6:8], extendSelection=False)

    e1 = p.edges
    p.MergeEdges(edgeList = e1[9:11], extendSelection=False)
    
    e = p.edges
    p.MergeEdges(edgeList = e[4:5]+e[8:9], extendSelection=False)
   
    e1 = p.edges
    p.MergeEdges(edgeList = e1[15:17], extendSelection=False)
    
    e = p.edges
    p.MergeEdges(edgeList = e[27:29], extendSelection=False)
    
    e1 = p.edges
    p.MergeEdges(edgeList = e1[25:27], extendSelection=False)
    mymodel.parts['Pressure_Vessel'].checkGeometry()

    ### Create Disturbance zone
    f, e, d1 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=d1[3], sketchUpEdge=e[3], 
        sketchPlaneSide=SIDE1, origin=(0.0, RV+LN, 
        H/2))
    s = mymodel.ConstrainedSketch(name='__profile__', 
        sheetSize=15236.66, gridSpacing=380.91, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    # Increase size of disturbance
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(RN+3/2*xsi_N, 0.0))
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#416 ]', ), )
    f1, e1, d2 = p.faces, p.edges, p.datums
    p.PartitionFaceBySketchThruAll(sketchPlane=d2[3], sketchUpEdge=e1[3], 
        faces=pickedFaces, sketchPlaneSide=SIDE1, sketch=s)
    s.unsetPrimaryObject()
    del mymodel.sketches['__profile__']

    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=RV+xsi_N)
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#2290 ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[28], faces=pickedFaces)

    ## Create Sets for Vessel Part
    session.viewports['Viewport: 1'].view.fitView()
    f = p.faces

    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#8000000 #640 ]', ), )
    p.Set(edges=edges, name='Edge_Vessel_Bottom')

    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#440 #c ]', ), )
    p.Set(edges=edges, name='Edge_Nozzle_Free')

    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#6105100a ]', ), )
    p.Set(edges=edges, name='Edge_Disturbance_Seed')
    
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#10022004 ]', ), )
    p.Set(edges=edges, name='Edge_Radial_Seed')

    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#5d6f0 ]', ), )
    p.Set(faces=faces, name='Surf_Vessel_Section')

    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#2290f ]', ), )
    p.Set(faces=faces, name='Surf_Nozzle_Section')

    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#416f0 ]', ), )
    p.Set(faces=faces, name='Surf_FrontFaces_Output')

    s = p.faces
    side2Faces = s.getSequenceFromMask(mask=('[#1250 ]', ), )
    p.Surface(side2Faces=side2Faces, name='SSurf_Vessel_Disturb')

    s = p.faces
    side2Faces = s.getSequenceFromMask(mask=('[#20109 ]', ), )
    p.Surface(side2Faces=side2Faces, name='SSurf_Nozzle_Disturb')

    s = p.faces 
    side1Faces = s.getSequenceFromMask(mask=('[#7ffff ]', ), )
    p.Surface(side1Faces=side1Faces, name='SSurf_Vessel_Pressure')

    s = p.edges
    side1Edges = s.getSequenceFromMask(mask=('[#0 #18a0 ]', ), )
    p.Surface(side1Edges=side1Edges, name='SEdge_Vessel_Top')

    return flatA

def Assembly(mymodel, thickness):
    """
    Assembles the pressure vessel. 
    Pressure Vessel part includes both main vessel and nozzle, 
    while reinformecement pad is in a separate part.
    Create instance for both and correctly position them.
    """
    t_ves, t_pad, t_noz = thickness

    t_off = (t_ves + t_pad)/2       # Offset from vessel is the average thickness of both components
    a = mymodel.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mymodel.parts['Pressure_Vessel']
    a.Instance(name='Pressure_Vessel-1', part=p, dependent=ON)

    p = mymodel.parts['Pad']
    a.Instance(name='Pad-1', part=p, dependent=ON)

    a.translate(instanceList=('Pad-1', ), vector=(0.0, t_off, 0.0))

    # Rotate the whole assembly so it matches with given local coordinate system
    a1 = mymodel.rootAssembly
    a1.rotate(instanceList=('Pressure_Vessel-1', 'Pad-1'), axisPoint=(0.0, 0.0, 
        0.0), axisDirection=(0.0, 0.0, 100.0), angle=-90.0)

    pass

def Property(mymodel, thickness):
    
    t_ves, t_pad, t_noz = thickness

    # Create Material
    mymodel.Material(name='X5CrNi18-10')
    mymodel.materials['X5CrNi18-10'].Elastic(table=((194000.0, 0.3), 
        ))
    
    # Create Sections
    mymodel.HomogeneousShellSection(name='Vessel_Section', 
        preIntegrate=OFF, material='X5CrNi18-10', thicknessType=UNIFORM, 
        thickness=t_ves, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
    mymodel.HomogeneousShellSection(name='Nozzle_Section', 
        preIntegrate=OFF, material='X5CrNi18-10', thicknessType=UNIFORM, 
        thickness=t_noz, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
    mymodel.HomogeneousShellSection(name='Pad_Section', 
        preIntegrate=OFF, material='X5CrNi18-10', thicknessType=UNIFORM, 
        thickness=t_pad, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
    
    # Assign Sections
    p = mymodel.parts['Pressure_Vessel']
    region = p.sets['Surf_Vessel_Section']
    p.SectionAssignment(region=region, sectionName='Vessel_Section', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    region = p.sets['Surf_Nozzle_Section']
    p.SectionAssignment(region=region, sectionName='Nozzle_Section', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    
    p = mymodel.parts['Pad']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    region = p.sets['Surf_Pad_Section']
    p.SectionAssignment(region=region, sectionName='Pad_Section', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    pass

def Step(mymodel):
    """
    Creates linear static loadstep and defines field output request
    """
    mymodel.StaticStep(name='Static_Vessel', previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        step='Static_Vessel')
    mymodel.fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', ), frequency=LAST_INCREMENT, position=AVERAGED_AT_NODES)
    
    # Uncomment if output is supposed to be for flat surface only
    regionDef=mdb.models['Model-1'].rootAssembly.allInstances['Pressure_Vessel-1'].sets['Surf_FrontFaces_Output']
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', ), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
        position=AVERAGED_AT_NODES)
    pass

def Interactions(mymodel):
    """
    Creates three interactions:
    TIE -> Reinforcement Pad to Nozzle
    TIE -> Reinforcement Pad to Vessel
    COUP -> Nozzle edge to center reference point
    """
    a = mymodel.rootAssembly
    a.regenerate()
    
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON)
    
    region1=a.instances['Pad-1'].sets['Edge_Pad_Outer']
    
    region2=a.instances['Pressure_Vessel-1'].surfaces['SSurf_Vessel_Disturb']
    mymodel.Tie(name='Pad_2_Vessel', main=region1, secondary=region2, 
        positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, 
        thickness=ON)
    
    region1=a.instances['Pad-1'].sets['Edge_Pad_Inner']
    
    region2=a.instances['Pressure_Vessel-1'].surfaces['SSurf_Nozzle_Disturb']
    mymodel.Tie(name='Pad_2_Nozzle', main=region1, secondary=region2, 
        positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, 
        thickness=ON)
    
    e1 = a.instances['Pressure_Vessel-1'].edges
    a.ReferencePoint(point=a.instances['Pressure_Vessel-1'].InterestingPoint(
        edge=e1[6], rule=CENTER))
    
    r1 = a.referencePoints
    refPoints1=(r1[6], )
    region1=a.Set(referencePoints=refPoints1, name='RP')
    
    region2=a.instances['Pressure_Vessel-1'].sets['Edge_Nozzle_Free']
    mymodel.Coupling(name='RBE3', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, 
        couplingType=DISTRIBUTING, 
        rotationalCouplingType=ROTATIONAL_STRUCTURAL, weightingMethod=UNIFORM, 
        localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)   

    pass

def Loads(mymodel, forces, moments, P, flatA, VR):
    FX, FY, FZ = forces
    MX, MY, MZ = moments
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, interactions=OFF, constraints=OFF, 
        engineeringFeatures=OFF)
    a = mymodel.rootAssembly
    region = a.sets['RP']

    # Nozzle Forces
    mymodel.ConcentratedForce(name='Nozzle_Forces', 
        createStepName='Static_Vessel', region=region, cf1=FX, 
        cf2=FY, cf3=FZ, distributionType=UNIFORM, field='', 
        localCsys=None)

    # Nozzle Moments
    mymodel.Moment(name='Nozzle_Moments', 
        createStepName='Static_Vessel', region=region, cm1=MX, cm2=MY, 
        cm3=MZ, distributionType=UNIFORM, field='', localCsys=None)
    
    # Top Edge Vertical Force
    mag = -P*flatA / (2*np.pi*VR)
    region = a.instances['Pressure_Vessel-1'].surfaces['SEdge_Vessel_Top']
    mdb.models['Model-1'].ShellEdgeLoad(name='Top_Edge_Force', 
        createStepName='Static_Vessel', region=region, magnitude=mag, 
        distributionType=UNIFORM, field='', localCsys=None, resultant=ON)

    # Internal Pressure
    region = a.instances['Pressure_Vessel-1'].surfaces['SSurf_Vessel_Pressure']
    mymodel.Pressure(name='Internal_Pressure', 
        createStepName='Static_Vessel', region=region, 
        distributionType=UNIFORM, field='', magnitude=P, amplitude=UNSET)
    
    # Clamp
    region = a.instances['Pressure_Vessel-1'].sets['Edge_Vessel_Bottom']
    mymodel.DisplacementBC(name='Bottom_Clamp', 
        createStepName='Static_Vessel', region=region, u1=0.0, u2=0.0, u3=0.0, 
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)

    pass

def Mesh(mymodel):
    elemType1 = mesh.ElemType(elemCode=S8R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=STRI65, elemLibrary=STANDARD)

    a = mymodel.rootAssembly
    a.regenerate()
    p = mymodel.parts['Pressure_Vessel']

    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#6105100a ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=12, constraint=FINER)

    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#10022004 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=30, constraint=FINER)

    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#63fff ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=QUAD, technique=STRUCTURED)
    p = mdb.models['Model-1'].parts['Pressure_Vessel']
    p.seedPart(size=150.0, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    faces = f.getSequenceFromMask(mask=('[#7ffff ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

    p = mymodel.parts['Pad']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#f ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=QUAD, technique=STRUCTURED)
    p.seedPart(size=30.0, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

    pass

def InputFile(mymodel, job_name='pvessel_000'):
    mdb.Job(name=job_name, model=mymodel.name, description='', type=ANALYSIS,
            memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, resultsFormat=ODB,
            numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    mdb.jobs[job_name].writeInput(consistencyChecking=OFF)

def PreProcess(VR, NR, H, LN, RP, thickness, forces, moments, P, jobname):
    Mdb()
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    mymodel = mdb.models[mdb.models.keys()[0]]

    VR, NR, H, LN, RP, thickness, xsi_N = validate_geom(VR, NR, H, LN, RP, thickness)
    flatA = Parts(mymodel, VR, NR, H, LN, RP, thickness)
    Assembly(mymodel, thickness)
    Property(mymodel, thickness)
    Step(mymodel)
    Interactions(mymodel)
    Loads(mymodel, forces, moments, P, flatA, VR)
    Mesh(mymodel)
    InputFile(mymodel, job_name=jobname)

    del mdb.models['Model-1']

    pass

# ======= Runing ======= #
def run_abaqus_job(folder_path, inp_file, job_name):
    try:
        print(f"Running {job_name} in {os.path.basename(folder_path)}...")
        run_job_script = os.path.join(os.path.dirname(__file__), "run_job.py")
        command = f'abaqus cae noGUI="{run_job_script}" -- "{inp_file}" "{job_name}"'
        subprocess.run(command, cwd=folder_path, shell=True, check=True)
        print(f"Done: {job_name}")
    except subprocess.CalledProcessError as e:
        print(f"Failed: {job_name} — {e}")

def parallel_exec(base_path, max_parallel_jobs):

    # Find all folders with .inp files and no .odb file
    jobs_to_run = []

    for folder_name in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder_name)
        if not os.path.isdir(folder_path):
            continue

        inp_files = [f for f in os.listdir(folder_path) if f.endswith('.inp')]
        if not inp_files:
            continue

        inp_file = inp_files[0]
        job_name = os.path.splitext(inp_file)[0]
        odb_path = os.path.join(folder_path, job_name + ".odb")

        if not os.path.exists(odb_path):
            jobs_to_run.append((folder_path, inp_file, job_name))
        else:
            print(f"{job_name}.odb already exists — skipping.")

    with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
        futures = [executor.submit(run_abaqus_job, folder, inp, job)
                   for folder, inp, job in jobs_to_run]
        for future in as_completed(futures):
            _ = future.result()

    print("All jobs finished.")


# ======= Post-Processing ======= #
def find_odb(parent_dir):
    odb_info = []
    # Loop through all entries in the parent directory
    for folder_name in os.listdir(parent_dir):
        folder_path = os.path.join(parent_dir, folder_name)
        # Check if it's a folder
        if os.path.isdir(folder_path):
            # Look for .odb files inside the folder
            for file_name in os.listdir(folder_path):
                if file_name.endswith('.odb'):
                    odb_path = os.path.join(folder_path, file_name)
                    odb_folder = os.path.dirname(odb_path)
                    odb_info.append((odb_path, odb_folder))
    return odb_info

def parse_inp(inp_path):
    from collections import defaultdict
    import re

    with open(inp_path, 'r') as f:
        lines = f.readlines()

    elset_to_thickness = {}
    elset_elements = defaultdict(list)
    elements_to_nodes = {}

    i = 0
    while i < len(lines):
        line = lines[i].strip().lower()

        # Parse *Shell Section
        if line.startswith("*shell section"):
            elset_match = re.search(r'elset=([\w\-]+)', line, re.IGNORECASE)
            if elset_match:
                elset_name = elset_match.group(1)
                thickness_line = lines[i + 1].strip()
                try:
                    thickness = float(thickness_line.split(',')[0])
                    elset_to_thickness[elset_name] = thickness
                except ValueError:
                    print(f"Warning: could not parse thickness on line {i+1}")
                i += 2
                continue

        # Parse *Elset (handle both normal and generate)
        elif line.startswith("*elset"):
            elset_match = re.search(r'elset=([\w\-]+)', line, re.IGNORECASE)
            generate = "generate" in line
            if elset_match:
                elset_name = elset_match.group(1)
                i += 1
                while i < len(lines) and not lines[i].strip().startswith("*"):
                    entries = [e.strip() for e in lines[i].strip().split(",") if e.strip()]
                    if generate and len(entries) == 3:
                        start, end, step = map(int, entries)
                        elset_elements[elset_name].extend(range(start, end + 1, step))
                    else:
                        elset_elements[elset_name].extend(map(int, entries))
                    i += 1
                continue

        # Parse *Element block (accept shell types like S4, S4R, S8R, etc.)
        elif line.startswith("*element") and "type=s" in line:
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = [int(p) for p in lines[i].strip().split(",") if p.strip()]
                if len(parts) >= 2:
                    elem_id = parts[0]
                    node_ids = parts[1:]
                    elements_to_nodes[elem_id] = node_ids
                i += 1
            continue

        i += 1  # Default increment

    # Now assign thickness per node
    node_thickness = defaultdict(list)
    for elset_name, element_ids in elset_elements.items():
        thickness = elset_to_thickness.get(elset_name)
        if thickness is None:
            print(f"Warning: No thickness found for elset '{elset_name}'")
            continue
        for elem_id in element_ids:
            node_ids = elements_to_nodes.get(elem_id, [])
            for node_id in node_ids:
                node_thickness[node_id].append(thickness)

    # Average if needed
    node_thickness_avg = {
        node: sum(t_list) / len(t_list)
        for node, t_list in node_thickness.items() if t_list
    }

    return node_thickness_avg

def PostProcess(odb_f, odb_dir, purge=0):
    """
    Purge==0 keeps everything
    Purge==1 keeps .inp, .obd and .csv
    Purge==2 keeps .csv only
    """
    odb = openOdb(path=odb_f)
    try:
        step_name = list(odb.steps.keys())[0]
        step = odb.steps[step_name]
        last_frame = step.frames[-1]

        stress = last_frame.fieldOutputs['S'].getSubset(position=NODAL)
        assembly = odb.rootAssembly

        base_filename = os.path.splitext(os.path.basename(odb_f))[0]
        csv_filename = os.path.join(odb_dir, f"{base_filename}.csv")
        inp_filename = os.path.join(odb_dir, f"{base_filename}.inp")

        # Parse .inp file for thickness
        if os.path.exists(inp_filename):
            thickness_map = parse_inp(inp_filename)
        else:
            print(f"Warning: .inp file not found: {inp_filename}")
            thickness_map = {}

        with open(csv_filename, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerow(['Part Name', 'Node ID', 'X', 'Y', 'Z', 'Thickness', 'S11', 'S22', 'VonMises'])

            written_nodes = set()

            for instance_name, instance in assembly.instances.items():
                coord_dict = {node.label: node.coordinates for node in instance.nodes}
                # Filter stress values by instance
                stress_subset = stress.getSubset(region=instance)

                for value in stress_subset.values:
                    node_label = value.nodeLabel
                    if (instance_name, node_label) not in written_nodes:
                        coords = coord_dict.get(node_label)
                        if coords is not None:  # Safe check
                            thickness = thickness_map.get(node_label, 0.0)
                            part_name = instance_name
                            node_id = node_label
                            
                            writer.writerow([
                                part_name,
                                node_id,
                                coords[0], coords[1], coords[2],
                                thickness,
                                value.data[0], value.data[1],
                                value.mises
                            ])
                            written_nodes.add((instance_name, node_label))
    except:
        pass

    odb.close()

    # Determine which files to keep
    if purge == 1:
        extensions_to_keep = {'.inp', '.odb', '.csv'}
    elif purge == 2:
        extensions_to_keep = {'.csv'}
        os.remove(odb_f)

    for file_name in os.listdir(odb_dir):
        file_path = os.path.join(odb_dir, file_name)
        if os.path.isfile(file_path):
            ext = os.path.splitext(file_name)[1].lower()
            if ext not in extensions_to_keep:
                try:
                    os.remove(file_path)
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Could not delete {file_path}: {e}")