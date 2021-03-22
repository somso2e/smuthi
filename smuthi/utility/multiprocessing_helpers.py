class ProcessesCluster:
    def __init__(self, max_processes):
        self.processes = []
        self.max_processes = max_processes


    def can_add_new_process(self):
        return len(self.processes) < self.max_processes


    def add_process(self, process):
        self.processes.append(process)


    def execute(self):
        for p in self.processes:
            p.start()

        for p in self.processes:
            p.join()


def distribute_processes_into_clusters(processes, max_processes_in_cluster):
    clusters = []
    cluster = ProcessesCluster(max_processes_in_cluster)
    clusters.append(cluster)

    for p in processes:
        if not cluster.can_add_new_process():
            cluster = ProcessesCluster(max_processes_in_cluster)
            clusters.append(cluster)

        cluster.add_process(p)

    return clusters