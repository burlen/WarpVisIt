import os
import sys
import visit_info
import parallel
import multiprocessing #Needed to start visit in a separate process
sys.path.append( visit_info.visitdir_sitepackages )
sys.path.append( visit_info.visitdir_bin )
#sys.path.append("/global/common/hopper2/graphics/visit/current/linux-x86_64/lib/site-packages")
#sys.path.append("/global/common/hopper2/graphics/visit/current/linux-x86_64/bin")

#Check if we are running within VisIt's CLI
isCLI = True
if 'VisItException' not in dir():
    if sys.version_info[1] < 5:
        print "\nRequires Python version 2.5 or greater\n"
        sys.exit()
    isCLI=False

#If we are in external python then import visit and launch
if not isCLI:
    from visit import *



def start_visit_viewer( simFile, vis_function,  visit_viewer_rank=0 , visit_viewer_rank_shared = True , visit_viewer_show_window =False , visit_viewer_debug=-1, visit_viewer_args=[]) :
        
    if visit_viewer_rank == parallel.get_rank() :
        print "Starting VisIt viewer"
        if visit_viewer_rank_shared :
            #We launch and run the viewer in a separate process so that it does not interfere with the simulation process"
            tViewer = multiprocessing.Process(target=run_visualization , args=(simFile, vis_function, visit_viewer_show_window, visit_viewer_debug, visit_viewer_args))
            tViewer.daemon = True #By making the viewer process a daemon allows us to call exit
            tViewer.start()
        else :
            run_visualization(simFile, vis_function, visit_viewer_show_window, visit_viewer_debug, visit_viewer_args )

def run_visualization(simFile , visFunction,  visit_viewer_show_window=False, visit_viewer_debug=-1, visit_viewer_args=[] ) :
    """Open the VisIt viewer, connect to the simulation and then execute the visualization defined by the visFunction. 
    """
    #Launch the viewer 
    launch_viewer(visit_viewer_show_window, visit_viewer_debug, visit_viewer_args  )
    #Connect to the simulation
    print "VisIt Viewer: open simulation"
    visit.OpenDatabase( simFile )
    print "VisIt Viewer: connected to the simulation"
    #Run the visualization script
    visFunction(simFile)


def launch_viewer( show_window=False , debug=-1 , viewer_args = [] ) :
    """Launch the visit viewer and add the dynamically loaded functions to the global scope of this module.
    
       Keyword arguments:
       show_window : Boolean indicating whether the VisIt viewer window should be shown.
                     Default value is False (add -nowin launch argument). If set to true
                     the -nowin launch argument will be omitted and ShowAllWindows() will
                     be called after the Launch is complete.
       debug : Integer indicating which debug level should be used for VisIt. Default is
               -1, disabling the debug output. Set to 1,2,3,4,5, to indicate which debug
               output should be provided.
       viewer_args : List of strings with additional input arguments for the VisIt 
                     viewer.
    """
    
    if not isCLI :
        #Add the launch arguments
        if not show_window :
            AddArgument( "-nowin" )
        if debug>0 :
            AddArgument( "-debug" )
            AddArgument( str(debug) )
        for a in viewer_args :
            AddArgument( a )
        
        #Launch the viewer
        Launch()
        #LaunchPySide()
    
        #VisIt has dynamically added a large range of functions
        #The following makes sure that we have all functions exposed globally
        import visit 
        for k, v in locals().items():
            globals()[k] = v
        
        #Make sure that all windows are shown if requested
        if show_window :
            visit.ShowAllWindows()
        print "Launched VisIt viewer"
        




