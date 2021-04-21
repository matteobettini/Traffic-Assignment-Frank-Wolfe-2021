import math
import time
import heapq
from scipy.optimize import fsolve

from network_import import *
from utils import PathUtils


class FlowTransportNetwork:
    tripSet = {}
    zoneSet = {}
    linkSet = {}
    nodeSet = {}
    originZones = {}

    @classmethod
    def reset(cls):
        cls.tripSet = {}
        cls.zoneSet = {}
        cls.linkSet = {}
        cls.nodeSet = {}
        cls.originZones = {}


class Zone:
    def __init__(self, zoneId: str):
        self.zoneId = zoneId

        self.lat = 0
        self.lon = 0
        self.destList = []  # list of zone ids (strs)


class Node:
    """
    This class has attributes associated with any node
    """

    def __init__(self, nodeId: str):
        self.Id = nodeId

        self.lat = 0
        self.lon = 0

        self.outLinks = []  # list of node ids (strs)
        self.inLinks = []  # list of node ids (strs)

        # For Dijkstra
        self.label = np.inf
        self.pred = None


class Link:
    """
    This class has attributes associated with any link
    """

    def __init__(self,
                 init_node: str,
                 term_node: str,
                 capacity: float,
                 length: float,
                 fft: float,
                 b: float,
                 power: float,
                 speed_limit: float,
                 toll: float,
                 linkType
                 ):
        self.init_node = init_node
        self.term_node = term_node
        self.capacity = float(capacity)  # veh per hour
        self.length = float(length)  # Length
        self.fft = float(fft)  # Free flow travel time (min)
        self.beta = float(power)
        self.alpha = float(b)
        self.speedLimit = float(speed_limit)
        self.toll = float(toll)
        self.linkType = linkType

        self.flow = 0.0
        self.cost = self.fft


class Demand:
    def __init__(self,
                 init_node: str,
                 term_node: str,
                 demand: float
                 ):
        self.fromZone = init_node
        self.toNode = term_node
        self.demand = float(demand)


def DijkstraHeap(origin):
    """
    Calcualtes shortest path from an origin to all other destinations.
    The labels and preds are stored in node instances.
    """
    for n in FlowTransportNetwork.nodeSet:
        FlowTransportNetwork.nodeSet[n].label = np.inf
        FlowTransportNetwork.nodeSet[n].pred = None
    FlowTransportNetwork.nodeSet[origin].label = 0.0
    FlowTransportNetwork.nodeSet[origin].pred = None
    SE = [(0, origin)]
    while SE:
        currentNode = heapq.heappop(SE)[1]
        currentLabel = FlowTransportNetwork.nodeSet[currentNode].label
        for toNode in FlowTransportNetwork.nodeSet[currentNode].outLinks:
            link = (currentNode, toNode)
            newNode = toNode
            newPred = currentNode
            existingLabel = FlowTransportNetwork.nodeSet[newNode].label
            newLabel = currentLabel + FlowTransportNetwork.linkSet[link].cost
            if newLabel < existingLabel:
                heapq.heappush(SE, (newLabel, newNode))
                FlowTransportNetwork.nodeSet[newNode].label = newLabel
                FlowTransportNetwork.nodeSet[newNode].pred = newPred


def BPRcostFunction(optimal: bool,
                    fft: float,
                    alpha: float,
                    flow: float,
                    capacity: float,
                    beta: float,
                    length: float,
                    maxSpeed: float
                    ) -> float:
    if optimal:
        return fft * (1 + (alpha * math.pow((flow * 1.0 / capacity), beta)) * (beta + 1))
    return fft * (1 + alpha * math.pow((flow * 1.0 / capacity), beta))


def constantCostFunction(optimal: bool,
                         fft: float,
                         alpha: float,
                         flow: float,
                         capacity: float,
                         beta: float,
                         length: float,
                         maxSpeed: float
                         ) -> float:
    if optimal:
        return fft + flow
    return fft


def greenshieldsCostFunction(optimal: bool,
                             fft: float,
                             alpha: float,
                             flow: float,
                             capacity: float,
                             beta: float,
                             length: float,
                             maxSpeed: float
                             ) -> float:
    if optimal:
        return (length * (capacity ** 2)) / (maxSpeed * (capacity - flow) ** 2)
    return length / (maxSpeed * (1 - (flow / capacity)))


def updateTravelTime(optimal: bool = False, costFunction=BPRcostFunction):
    """
    This method updates the travel time on the links with the current flow
    """
    for l in FlowTransportNetwork.linkSet:
        FlowTransportNetwork.linkSet[l].cost = costFunction(optimal,
                                                            FlowTransportNetwork.linkSet[l].fft,
                                                            FlowTransportNetwork.linkSet[l].alpha,
                                                            FlowTransportNetwork.linkSet[l].flow,
                                                            FlowTransportNetwork.linkSet[l].capacity,
                                                            FlowTransportNetwork.linkSet[l].beta,
                                                            FlowTransportNetwork.linkSet[l].length,
                                                            FlowTransportNetwork.linkSet[l].speedLimit
                                                            )


def findAlpha(x_bar, optimal: bool = False, costFunction=BPRcostFunction):
    """
    This uses unconstrained optimization to calculate the optimal step size required
    for Frank-Wolfe Algorithm
    """

    def df(alpha):
        sum_derivative = 0  # this line is the derivative of the objective function.
        for l in FlowTransportNetwork.linkSet:
            tmpFlow = alpha * x_bar[l] + (1 - alpha) * FlowTransportNetwork.linkSet[l].flow
            tmpCost = costFunction(optimal,
                                   FlowTransportNetwork.linkSet[l].fft,
                                   FlowTransportNetwork.linkSet[l].alpha,
                                   tmpFlow,
                                   FlowTransportNetwork.linkSet[l].capacity,
                                   FlowTransportNetwork.linkSet[l].beta,
                                   FlowTransportNetwork.linkSet[l].length,
                                   FlowTransportNetwork.linkSet[l].speedLimit
                                   )
            sum_derivative = sum_derivative + (x_bar[l] - FlowTransportNetwork.linkSet[l].flow) * tmpCost
        return sum_derivative

    sol = fsolve(df, np.array([0.5]))
    return max(0, min(1, sol[0]))


def tracePreds(dest):
    """
    This method traverses predecessor nodes in order to create a shortest path
    """
    prevNode = FlowTransportNetwork.nodeSet[dest].pred
    spLinks = []
    while prevNode is not None:
        spLinks.append((prevNode, dest))
        dest = prevNode
        prevNode = FlowTransportNetwork.nodeSet[dest].pred
    return spLinks


def loadAON(computeXbar: bool = True):
    """
    This method produces auxiliary flows for all or nothing loading.
    """
    x_bar = {l: 0.0 for l in FlowTransportNetwork.linkSet}
    SPTT = 0.0
    for r in FlowTransportNetwork.originZones:
        DijkstraHeap(r)
        for s in FlowTransportNetwork.zoneSet[r].destList:
            dem = FlowTransportNetwork.tripSet[r, s].demand

            SPTT = SPTT + FlowTransportNetwork.nodeSet[s].label * dem

            if computeXbar and r != s:
                for spLink in tracePreds(s):
                    x_bar[spLink] = x_bar[spLink] + dem

    return SPTT, x_bar

def readDemand(demand_df: pd.DataFrame):
    for index, row in demand_df.iterrows():

        init_node = str(int(row["init_node"]))
        term_node = str(int(row["term_node"]))
        demand = row["demand"]
        if demand <= 0:
            continue

        FlowTransportNetwork.tripSet[init_node, term_node] = Demand(init_node, term_node, demand)
        if init_node not in FlowTransportNetwork.zoneSet:
            FlowTransportNetwork.zoneSet[init_node] = Zone(init_node)
        if term_node not in FlowTransportNetwork.zoneSet:
            FlowTransportNetwork.zoneSet[term_node] = Zone(term_node)
        if term_node not in FlowTransportNetwork.zoneSet[init_node].destList:
            FlowTransportNetwork.zoneSet[init_node].destList.append(term_node)

    print(len(FlowTransportNetwork.tripSet), "OD pairs")
    print(len(FlowTransportNetwork.zoneSet), "OD zones")


def readNetwork(network_df: pd.DataFrame):
    for index, row in network_df.iterrows():

        init_node = str(int(row["init_node"]))
        term_node = str(int(row["term_node"]))
        capacity = row["capacity"]
        length = row["length"]
        free_flow_time = row["free_flow_time"]
        b = row["b"]
        power = row["power"]
        speed = row["speed"]
        toll = row["toll"]
        link_type = row["link_type"]

        FlowTransportNetwork.linkSet[init_node, term_node] = Link(init_node=init_node,
                                                                  term_node=term_node,
                                                                  capacity=capacity,
                                                                  length=length,
                                                                  fft=free_flow_time,
                                                                  b=b,
                                                                  power=power,
                                                                  speed_limit=speed,
                                                                  toll=toll,
                                                                  linkType=link_type
                                                                  )
        if init_node not in FlowTransportNetwork.nodeSet:
            FlowTransportNetwork.nodeSet[init_node] = Node(init_node)
        if term_node not in FlowTransportNetwork.nodeSet:
            FlowTransportNetwork.nodeSet[term_node] = Node(term_node)
        if term_node not in FlowTransportNetwork.nodeSet[init_node].outLinks:
            FlowTransportNetwork.nodeSet[init_node].outLinks.append(term_node)
        if init_node not in FlowTransportNetwork.nodeSet[term_node].inLinks:
            FlowTransportNetwork.nodeSet[term_node].inLinks.append(init_node)

    print(len(FlowTransportNetwork.nodeSet), "nodes")
    print(len(FlowTransportNetwork.linkSet), "links")



def assignment_loop(algorithm: str = "FW",
                    systemOptimal: bool = False,
                    costFunction=BPRcostFunction,
                    accuracy: float = 0.001,
                    maxIter: int = 1000,
                    maxTime: int = 60):
    """
    For explaination of the algorithm see Chapter 7 of:
    https://sboyles.github.io/blubook.html
    PDF:
    https://sboyles.github.io/teaching/ce392c/book.pdf
    """

    iteration_number = 1
    gap = np.inf
    assignmentStartTime = time.time()

    # Check if desired accuracy is reached
    while gap > accuracy:

        # Get x_bar throug all-or-nothing assignment
        _, x_bar = loadAON()

        if algorithm == "MSA" or iteration_number == 1:
            alpha = (1 / iteration_number)
        elif algorithm == "FW":
            # If using Frank-Wolfe determine the step size alpha by solving a nonlinear equation
            alpha = findAlpha(x_bar,
                              optimal=systemOptimal,
                              costFunction=costFunction)
        else:
            print("Terminating the program.....")
            print("The solution algorithm ", algorithm, " does not exist!")
            return

        # Apply flow improvement
        for l in FlowTransportNetwork.linkSet:
            FlowTransportNetwork.linkSet[l].flow = alpha * x_bar[l] + (1 - alpha) * FlowTransportNetwork.linkSet[l].flow

        # Compute the new travel time
        updateTravelTime(optimal=systemOptimal,
                         costFunction=costFunction)

        # Compute the relative gap
        SPTT, _ = loadAON(computeXbar=False)
        TSTT = sum([FlowTransportNetwork.linkSet[a].flow * FlowTransportNetwork.linkSet[a].cost for a in
                    FlowTransportNetwork.linkSet])
        gap = (TSTT / SPTT) - 1
        assert gap >= 0

        iteration_number += 1
        if iteration_number > maxIter:
            print("The assignment did not converge to the desired gap and the max number of iterations has been reached")
            print("Assignment took", round(time.time() - assignmentStartTime, 5), "seconds")
            print("Current gap:", round(gap, 5))
            return
        if time.time() - assignmentStartTime > maxTime:
            print("The assignment did not converge to the desired gap and the max time limit has been reached")
            print("Assignment did ", iteration_number, "iterations")
            print("Current gap:", round(gap, 5))
            return

    print("Assignment converged in ", iteration_number, "iterations")
    print("Assignment took", round(time.time() - assignmentStartTime, 5), "seconds")
    print("Current gap:", round(gap, 5))


def writeResults(output_file, costFunction=BPRcostFunction, systemOptimal: bool = False):
    outFile = open(output_file, "w")
    TSTT = round(sum([FlowTransportNetwork.linkSet[a].flow * costFunction(False,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].fft,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].alpha,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].flow,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].capacity,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].beta,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].length,
                                                                          FlowTransportNetwork.linkSet[
                                                                              a].speedLimit
                                                                          ) for a in
                      FlowTransportNetwork.linkSet]), 3)
    print("\nTotal system travel time:", f'{TSTT} secs')
    tmpOut = "Total Travel Time:\t" + str(TSTT)
    outFile.write(tmpOut + "\n")
    tmpOut = "Cost function used:\t" + BPRcostFunction.__name__
    outFile.write(tmpOut + "\n")
    tmpOut = ["User equilibrium (UE) or system optimal (SO):\t"] + ["SO" if systemOptimal else "UE"]
    outFile.write("".join(tmpOut) + "\n\n")
    tmpOut = "init_node\tterm_node\tflow\ttravelTime"
    outFile.write(tmpOut + "\n")
    for i in FlowTransportNetwork.linkSet:
        tmpOut = str(FlowTransportNetwork.linkSet[i].init_node) + "\t" + str(
            FlowTransportNetwork.linkSet[i].term_node) + "\t" + str(
            FlowTransportNetwork.linkSet[i].flow) + "\t" + str(costFunction(False,
                                                                            FlowTransportNetwork.linkSet[i].fft,
                                                                            FlowTransportNetwork.linkSet[i].alpha,
                                                                            FlowTransportNetwork.linkSet[i].flow,
                                                                            FlowTransportNetwork.linkSet[i].capacity,
                                                                            FlowTransportNetwork.linkSet[i].beta,
                                                                            FlowTransportNetwork.linkSet[i].length,
                                                                            FlowTransportNetwork.linkSet[i].speedLimit
                                                                            ))
        outFile.write(tmpOut + "\n")
    outFile.close()


def computeAssingment(net_file: str,
                      demand_file: str = None,
                      algorithm: str = "FW",  # FW or MSA
                      costFunction=BPRcostFunction,
                      systemOptimal: bool = False,
                      accuracy: float = 0.0001,
                      maxIter: int = 1000,
                      maxTime: int = 60,
                      results_file: str = None
                      ) -> None:
    """
    This is the main function to compute the user equilibrium UE (default) or system optimal (SO) traffic assignment
    All the networks present on https://github.com/bstabler/TransportationNetworks following the tntp format can be loaded

    :param net_file: Name of the network (net) file following the tntp format (see https://github.com/bstabler/TransportationNetworks)
    :param demand_file: Name of the demand (trips) file following the tntp format (see https://github.com/bstabler/TransportationNetworks)
    :param algorithm:
           - "FW": Frank-Wolfe algorithm (see https://en.wikipedia.org/wiki/Frank%E2%80%93Wolfe_algorithm)
           - "MSA": Method of successive averages
           For more information on how the algorithms work see https://sboyles.github.io/teaching/ce392c/book.pdf
    :param costFunction: Which cost function to use to compute travel time on edges, currently available functions are:
           - BPRcostFunction (see https://rdrr.io/rforge/travelr/man/bpr.function.html)
           - greenshieldsCostFunction (see Greenshields, B. D., et al. "A study of traffic capacity." Highway research board proceedings. Vol. 1935. National Research Council (USA), Highway Research Board, 1935.)
           - constantCostFunction
    :param systemOptimal: Wheather to compute the system optimal flows instead of the user equilibrium
    :param accuracy: Desired assignment precision gap
    :param maxIter: Maximum nuber of algorithm iterations
    :param maxTime: Maximum seconds allowed for the assignment
    :param results_file: Name of the desired file to write the results,
           by default the result file is saved with the same name as the input network with the suffix "_flow.tntp" in the same folder
    :return: None
    """

    if demand_file is None:
        demand_file = '_'.join(net_file.split("_")[:-1] + ["trips.tntp"])

    readStart = time.time()

    net_name = net_file.split("/")[-1].split("_")[0]
    print(f"Loading network {net_name}...")

    net_df, demand_df = import_network(
        net_file,
        demand_file)

    readDemand(demand_df)
    readNetwork(net_df)

    FlowTransportNetwork.originZones = set([k[0] for k in FlowTransportNetwork.tripSet])

    print("Network", net_name, "loaded")
    print("Reading the network data took", round(time.time() - readStart, 2), "secs\n")
    print("Computing assignment...")
    assignment_loop(algorithm=algorithm, systemOptimal=systemOptimal, costFunction=costFunction, accuracy=accuracy, maxIter=maxIter, maxTime=maxTime)

    if results_file is None:
        results_file = '_'.join(net_file.split("_")[:-1] + ["flow.tntp"])

    writeResults(results_file,
                 costFunction=costFunction,
                 systemOptimal=systemOptimal)

    FlowTransportNetwork.reset()


if __name__ == '__main__':
    computeAssingment(net_file=str(PathUtils.sioux_falls_net_file),
                      algorithm="FW",
                      costFunction=BPRcostFunction,
                      systemOptimal=False,
                      accuracy=0.0001,
                      maxIter=1000,
                      maxTime=60)
