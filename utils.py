from pathlib import Path


def get_project_root() -> Path:
    return Path(__file__).resolve().parent.absolute()


class PathUtils:
    """
    Utils for file names
    """

    # FOLDERS

    # Input netwroks
    input_networks_folder = get_project_root() / "InputNetworks"

    # FILES

    # Network files
    anaheim_net_file = input_networks_folder / "Anaheim_net.tntp"
    barcelona_net_file = input_networks_folder / "Barcelona_net.tntp"
    braess_net_file = input_networks_folder / "Braess_net.tntp"
    chicago_net_file = input_networks_folder / "ChicagoSketch_net.tntp"
    sioux_falls_net_file = input_networks_folder / "SiouxFalls_net.tntp"
    winnipeg_net_file = input_networks_folder / "Winnipeg_net.tntp"

    # Demand files
    # anaheim_demand_file = input_networks_folder / "Anaheim_trips.tntp"
    # barcelona_demand_file = input_networks_folder / "Barcelona_trips.tntp"
    # braess_demand_file = input_networks_folder / "Braess_trips.tntp"
    # chicago_demand_file = input_networks_folder / "ChicagoSketch_trips.tntp"
    # sioux_falls_demand_file = input_networks_folder / "SiouxFalls_trips.tntp"
    # winnipeg_demand_file = input_networks_folder / "Winnipeg_trips.tntp"


