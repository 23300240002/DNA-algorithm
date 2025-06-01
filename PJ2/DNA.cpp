#include <cstdio>
#include <vector>
#include <deque>
#include <unordered_map>
#include <unordered_set>

// 匹配状态节点，type表示状态类型，i为reference位置，j为query位置
struct NodeEntry {
    char type;
    unsigned short i, j;
    bool operator==(const NodeEntry& other) const {
        return type == other.type && i == other.i && j == other.j;
    }
};

// 自定义哈希函数，用于unordered_map和unordered_set
namespace std {
template <>
struct hash<NodeEntry> {
    size_t operator()(const NodeEntry& entry) const {
        return (size_t(entry.i) << 16) | entry.j;
    }
};
}

FILE* ref_file;
FILE* query_file;
FILE* plot_file;
FILE* dbg_file;

// 距离信息，dis为当前距离，entry为前驱节点
struct dis_t {
    int dis;
    NodeEntry entry;
};

std::vector<char> reference, reference_r, query;
std::unordered_map<NodeEntry, dis_t, std::hash<NodeEntry>> dis;
std::unordered_set<NodeEntry, std::hash<NodeEntry>> vis;
std::deque<NodeEntry> queue;

// 互补碱基转换
char tr(char c) {
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

int main(int argc, char* argv[]) {
    // 读取输入文件
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <reference_file> <query_file> \n", argv[0]);
        return 1;
    }
    ref_file = fopen(argv[1], "r");
    if (ref_file == NULL) {
        fprintf(stderr, "Error opening reference file\n");
        return 1;
    }
    query_file = fopen(argv[2], "r");
    if (query_file == NULL) {
        fprintf(stderr, "Error opening query file\n");
        fclose(ref_file);
        return 1;
    }
    plot_file = fopen("plot.txt", "w");
    if (plot_file == NULL) {
        fprintf(stderr, "Error opening plot file\n");
        fclose(ref_file);
        fclose(query_file);
        return 1;
    }
    dbg_file = fopen("debug.txt", "w");
    if (dbg_file == NULL) {
        fprintf(stderr, "Error opening debug file\n");
        fclose(ref_file);
        fclose(query_file);
        fclose(plot_file);
        return 1;
    }

    reference.clear();
    reference_r.clear();
    query.clear();
    // 占位符，方便下标从1开始
    reference.push_back('X');
    query.push_back('X');
    reference_r.push_back('Y');
    char c;
    // 读取reference序列，同时生成其互补序列
    while (true) {
        c = fgetc(ref_file);
        if (c == EOF) break;
        if (c == '\n') break;
        reference.push_back(c);
        reference_r.push_back(tr(c));
    }
    // 读取query序列
    while (true) {
        c = fgetc(query_file);
        if (c == EOF) break;
        if (c == '\n') break;
        query.push_back(c);
    }

    // 初始化起点，type为'C'，表示可以从任意reference位置开始匹配
    dis[NodeEntry{'C', 0, 0}] = {0, NodeEntry{'C', 0, 0}};
    queue.push_back(NodeEntry{'C', 0, 0});
    unsigned short tgt_j = query.size() - 1;

    // 广度优先搜索，寻找最优匹配路径
    while (!queue.empty()) {
        NodeEntry entry = queue.front();
        // 如果已经匹配到query末尾，结束
        if (entry.j == tgt_j) {
            printf("Found a path with distance %d\n", dis[entry].dis);
            break;
        }
        queue.pop_front();
        if (vis.find(entry) != vis.end()) {
            continue;
        }
        vis.insert(entry);
        int distance = dis[entry].dis;

        // C状态：可以选择从reference正向或反向互补任意位置开始新的匹配
        if (entry.type == 'C') {
            for (int i = 0; i <= (int)reference.size(); i++) {
                NodeEntry next_entry{'A', (unsigned short)i, entry.j};
                // C++17没有contains，改用find
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                    dis[next_entry] = {distance + 1, entry};
                    queue.push_back(next_entry);
                }
                NodeEntry next_entry2{'B', (unsigned short)i, entry.j};
                if (dis.find(next_entry2) == dis.end() || dis[next_entry2].dis > distance + 1) {
                    dis[next_entry2] = {distance + 1, entry};
                    queue.push_back(next_entry2);
                }
            }
            printf("Reach query %d, distance: %d\n", entry.j, distance);
        }
        // A状态：reference正向与query比对
        else if (entry.type == 'A') {
            unsigned short i = entry.i;
            unsigned short j = entry.j;
            // 尝试匹配下一个碱基
            if (i + 1 <= reference.size() - 1 && j + 1 <= query.size() - 1) {
                NodeEntry next_entry{'A', (unsigned short)(i + 1), (unsigned short)(j + 1)};
                int new_distance = (reference[i + 1] == query[j + 1]) ? 0 : 1; // 碱基相同则距离不变，否则+1
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + new_distance) {
                    dis[next_entry] = {distance + new_distance, entry};
                    if (new_distance == 0) {
                        queue.push_front(next_entry); // 优先扩展无错配的路径
                    } else {
                        queue.push_back(next_entry);
                    }
                }
            }
            // 跳过reference一个碱基（相当于删除）
            if (i + 1 <= reference.size() - 1) {
                NodeEntry next_entry{'A', (unsigned short)(i + 1), j};
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                    dis[next_entry] = {distance + 1, entry};
                    queue.push_back(next_entry);
                }
            }
            // 跳过query一个碱基（相当于插入）
            if (j + 1 <= query.size() - 1) {
                NodeEntry next_entry{'A', i, (unsigned short)(j + 1)};
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                    dis[next_entry] = {distance + 1, entry};
                    queue.push_back(next_entry);
                }
            }
            // 可以随时切换到C状态，重新选择锚点
            NodeEntry next_entry{'C', 0, j};
            if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                dis[next_entry] = {distance + 1, entry};
                queue.push_back(next_entry);
            }
        }
        // B状态：reference反向互补与query比对
        else if (entry.type == 'B') {
            unsigned short i = entry.i;
            unsigned short j = entry.j;
            // 反向互补匹配下一个碱基
            if (i >= 1 && j + 1 <= query.size() - 1) {
                NodeEntry next_entry{'B', (unsigned short)(i - 1), (unsigned short)(j + 1)};
                int new_distance = (reference_r[i - 1] == query[j + 1]) ? 0 : 1;
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + new_distance) {
                    dis[next_entry] = {distance + new_distance, entry};
                    if (new_distance == 0) {
                        queue.push_front(next_entry);
                    } else {
                        queue.push_back(next_entry);
                    }
                }
            }
            // 跳过reference一个碱基（反向）
            if (i >= 1) {
                NodeEntry next_entry{'B', (unsigned short)(i - 1), j};
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                    dis[next_entry] = {distance + 1, entry};
                    queue.push_back(next_entry);
                }
            }
            // 跳过query一个碱基
            if (j + 1 <= query.size() - 1) {
                NodeEntry next_entry{'B', i, (unsigned short)(j + 1)};
                if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                    dis[next_entry] = {distance + 1, entry};
                    queue.push_back(next_entry);
                }
            }
            // 可以随时切换到C状态，重新选择锚点
            NodeEntry next_entry{'C', 0, j};
            if (dis.find(next_entry) == dis.end() || dis[next_entry].dis > distance + 1) {
                dis[next_entry] = {distance + 1, entry};
                queue.push_back(next_entry);
            }
        }
    }

    // 回溯输出匹配路径，写入plot.txt
    NodeEntry entry = queue.front();
    int seg_end_i = entry.i, seg_end_j = entry.j;
    while (entry.type != 'C' || entry.i != 0 || entry.j != 0) {
        // 正向匹配段输出
        if (entry.type == 'A' && entry.type != dis[entry].entry.type) {
            fprintf(plot_file, "(%d,%d,%d,%d),\n", entry.j, seg_end_j, entry.i, seg_end_i);
            fflush(plot_file);
        }
        // 反向互补匹配段输出
        if (entry.type == 'B' && entry.type != dis[entry].entry.type) {
            fprintf(plot_file, "(%d,%d,%d,%d),\n", entry.j, seg_end_j, seg_end_i - 1, entry.i - 1);
            fflush(plot_file);
        }
        // C状态，更新当前段终点
        if (entry.type == 'C' && entry.type != dis[entry].entry.type) {
            seg_end_i = dis[entry].entry.i;
            seg_end_j = dis[entry].entry.j;
        }
        entry = dis[entry].entry;
    }
}