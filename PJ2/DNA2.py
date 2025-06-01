from collections import deque, namedtuple
from typing import List, Tuple, Dict, Optional
import edlib # 用于DNA序列比对的库，高效的编辑距离计算

# NodeEntry是一个命名元组，表示BFS中的节点状态，一个节点包含三部分信息：
# type: A表示查询序列，B表示反向互补序列，C表示参考序列
# i: 参考序列位置，j: 查询序列位置
# 例如，NodeEntry('A', 5, 7)表示：在query中位置7，reference中位置5，此时处于A状态
Node = namedtuple('Node', ['type', 'i', 'j'])

# 互补碱基转换
def transition(c):
    c = c.upper()
    transitions = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    if c not in transitions:
        raise ValueError(f"DNA序列不合法，出现: {c}")
    return transitions[c]

# 输入一个DNA序列，返回其反向互补序列
def rc(s):
    return ''.join(transition(c) for c in s)[::-1]

# 预处理输入序列，生成反向互补
def preprocess_sequences(reference, query):
    reference = 'P' + reference # 占位符保证真实数组从1开始
    query = 'P' + query
    reverse = 'Q' + ''.join(transition(c) for c in reference[1:]) # ''表示用空字符串连接
    n = len(reference) - 1
    m = len(query) - 1
    return reference, query, reverse, n, m

# 核心算法：根据当前状态扩展所有可能的下一个状态，并入队
def expand_state(node, reference, query, reverse, n, m, state_map, state_queue):
    distance = state_map[node][0] # 当前状态的最小步数
    if node.type == 'C': # 当前需要重新选择一个reference上的位置作为新的比对起点
        for i in range(n + 1): # reference一步步往前
            next_node = Node('A', i, node.j) # 尝试正向比对
            if next_node not in state_map or state_map[next_node][0] > distance + 1:
                # 如果新状态不在字典中，或者该状态的距离可以更小，则更新
                state_map[next_node] = (distance + 1, node)
                state_queue.append(next_node)
            next_node2 = Node('B', i, node.j) # 尝试反向比对
            if next_node2 not in state_map or state_map[next_node2][0] > distance + 1:
                state_map[next_node2] = (distance + 1, node)
                state_queue.append(next_node2)
    elif node.type == 'A': # 当前在正向比对状态
        i, j = node.i, node.j
        if i + 1 <= n and j + 1 <= m: # 参考序列和查询序列都往前走一步
            next_node = Node('A', i + 1, j + 1) # 尝试正向比对
            # 如果当前参考序列和查询序列的碱基相同，则不增加距离，否则增加1
            new_distance = 0 if reference[i + 1] == query[j + 1] else 1
            if next_node not in state_map or state_map[next_node][0] > distance + new_distance:
                state_map[next_node] = (distance + new_distance, node)
                if new_distance == 0: # 如果没有增加距离，必须优先入队，这样才能保证后续BFS时前面的状态已处理好！
                    state_queue.appendleft(next_node)
                else:
                    state_queue.append(next_node)
        if i + 1 <= n: # 尝试不走查询序列，只走参考序列，即发生一次缺失
            next_node = Node('A', i + 1, j)
            if next_node not in state_map or state_map[next_node][0] > distance + 1:
                state_map[next_node] = (distance + 1, node)
                state_queue.append(next_node)
        if j + 1 <= m: # 尝试不走参考序列，只走查询序列，即发生一次插入
            next_node = Node('A', i, j + 1)
            if next_node not in state_map or state_map[next_node][0] > distance + 1:
                state_map[next_node] = (distance + 1, node)
                state_queue.append(next_node)
        next_node = Node('C', 0, j) # 还可以从query当前位置j和reference的任意新位置（从0开始）重新开始比对
        if next_node not in state_map or state_map[next_node][0] > distance + 1:
            state_map[next_node] = (distance + 1, node)
            state_queue.append(next_node)
    elif node.type == 'B': # 当前在反向比对状态
        i, j = node.i, node.j
        if i >= 1 and j + 1 <= m: # 参考序列和查询序列都往前走一步
            next_node = Node('B', i - 1, j + 1) # 尝试反向比对
            new_distance = 0 if reverse[i - 1] == query[j + 1] else 1
            if next_node not in state_map or state_map[next_node][0] > distance + new_distance:
                state_map[next_node] = (distance + new_distance, node)
                if new_distance == 0:
                    state_queue.appendleft(next_node)
                else:
                    state_queue.append(next_node)
        if i >= 1: # 尝试不走查询序列，只走参考序列，即发生一次缺失
            next_node = Node('B', i - 1, j)
            if next_node not in state_map or state_map[next_node][0] > distance + 1:
                state_map[next_node] = (distance + 1, node)
                state_queue.append(next_node)
        if j + 1 <= m: # 尝试不走参考序列，只走查询序列，即发生一次插入
            next_node = Node('B', i, j + 1)
            if next_node not in state_map or state_map[next_node][0] > distance + 1:
                state_map[next_node] = (distance + 1, node)
                state_queue.append(next_node)
        next_node = Node('C', 0, j) # 还可以从query当前位置j和reference的任意新位置（从0开始）重新开始比对
        if next_node not in state_map or state_map[next_node][0] > distance + 1:
            state_map[next_node] = (distance + 1, node)
            state_queue.append(next_node)

# 主BFS循环，负责图遍历和距离记录
def bfs_search(reference, query, reverse, n, m):
    state_map = {}  # key: Node, value: (距离, 前驱节点)
    visited = set()
    state_queue = deque()
    start = Node('C', 0, 0)
    state_map[start] = (0, None)
    state_queue.append(start)

    # 新增：每个query位置的最优距离
    best_distance_at_j = dict()
    # 新增：全局最优终点距离
    min_distance_to_end = None
    # 剪枝阈值
    PRUNE_MARGIN = 2

    max_percent = -1
    max_query_j = -1
    found_entry = None
    end_nodes = []  # 新增

    while state_queue:
        node = state_queue.popleft()
        distance = state_map[node][0]

        # 剪枝1：只保留每个query位置的最优距离
        if node.j in best_distance_at_j:
            if distance > best_distance_at_j[node.j] + PRUNE_MARGIN:
                continue
        # 更新最优
        if node.j not in best_distance_at_j or distance < best_distance_at_j[node.j]:
            best_distance_at_j[node.j] = distance

        # 剪枝2：如果已知终点最优距离，丢弃明显更差的路径
        if min_distance_to_end is not None and distance > min_distance_to_end + PRUNE_MARGIN:
            continue

        # 进度输出
        if node.j > max_query_j:
            max_query_j = node.j
            percent = int(node.j * 100 / m)
            if percent > max_percent:
                max_percent = percent
                print(f"当前进展至 {percent}%（query位置 {node.j}/{m}）")
        if node.j == m:
            end_nodes.append(node)  # 修复：加入终点节点
            # found_entry = node
            # 记录最优终点距离
            if min_distance_to_end is None or distance < min_distance_to_end:
                min_distance_to_end = distance
            continue  # 不break，继续找更优的终点

        if node in visited:
            continue
        visited.add(node)
        expand_state(node, reference, query, reverse, n, m, state_map, state_queue)

    # 修复：选取所有j==m节点中距离最小的
    if end_nodes:
        found_entry = min(end_nodes, key=lambda n: state_map[n][0])
    else:
        found_entry = None
    return state_map, found_entry

# 回溯最优路径，生成最终结果
def backtrack_path(state_map, found_entry):
    start = Node('C', 0, 0) # 起始节点
    result = [] # 存储最后的所有路径
    if found_entry is None:
        print("未找到匹配路径")
        return result
    node = found_entry # 从找到的节点状态（即路径的终点）开始回溯，不断更新为其前驱
    seg_end_i, seg_end_j = node.i, node.j # 记录当前段的结束位置
    # 反向回溯，直到回到起始节点
    while node != start:
        prev_entry = state_map[node][1]
        if node.type == 'A' and (prev_entry is None or node.type != prev_entry.type):
            # 只保留长度>=30的段
            if seg_end_j - node.j >= 30:
                result.append((node.j, seg_end_j, node.i, seg_end_i))
        if node.type == 'B' and (prev_entry is None or node.type != prev_entry.type):
            if seg_end_j - node.j >= 30:
                result.append((node.j, seg_end_j, seg_end_i - 1, node.i - 1))
        if node.type == 'C' and (prev_entry is None or node.type != prev_entry.type):
            seg_end_i = prev_entry.i
            seg_end_j = prev_entry.j
        node = prev_entry
    result.reverse() # 保证输出顺序从前到后
    print("比对完成！")
    return result

# 计算编辑距离（来自ipynb的测试方法）
def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    A = ref[ref_st: ref_en]
    a = query[query_st: query_en]
    _a = rc(query[query_st: query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

# 筛选失配比例低的
def filter_segments(segments, ref, query, max_mismatch_ratio=0.1):
    filtered = []
    for seg in segments:
        query_st, query_en, ref_st, ref_en = seg
        seg_len = query_en - query_st
        if seg_len < 30:
            continue
        edit_dist = calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
        if seg_len == 0:
            continue
        mismatch_ratio = edit_dist / seg_len
        if mismatch_ratio <= max_mismatch_ratio: # 只保留性价比高片段
            filtered.append(seg)
    return filtered

# 主函数，分为三个步骤
def dna_align(reference, query): 
    # 1. 预处理输入序列，添加占位符并生成互补翻转序列
    reference, query, reverse, n, m = preprocess_sequences(reference, query)
    # 2. BFS搜索最短路径，记录每个状态的距离和前驱节点
    state_map, found_entry = bfs_search(reference, query, reverse, n, m)
    # 3. 回溯路径，生成最终结果
    result = backtrack_path(state_map, found_entry)
    # 4. 最后筛选“性价比高”的片段
    filtered = filter_segments(result, reference[1:], query[1:]) # 去掉占位符
    return filtered

if __name__ == "__main__":
    reference_file = "reference2.txt"
    query_file = "query2.txt"
    output_file = "answer2.txt"
    
    # 读取参考序列和查询序列
    with open(reference_file, "r") as f:
        reference = f.read().strip()
    with open(query_file, "r") as f:
        query = f.read().strip()
    
    # 调用比对函数并输出结果
    result = dna_align(reference, query)
    with open(output_file, "w") as f:
        for item in result:
            f.write(f"{item}\n")
    print(f"比对结果已写入 {output_file}，请助教查阅。")