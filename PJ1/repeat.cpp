#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
using namespace std;

// 计算字符串的互补翻转序列
string get_reverse_complement(const string& s) {
    string res = s;
    reverse(res.begin(), res.end());
    for (char& c : res) {
        if (c == 'A') c = 'T';
        else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G';
        else if (c == 'G') c = 'C';
    }
    return res;
}

int main() {
    string reference, query;
    cout << "Please input the reference sequence:" << endl;
    cin >> reference;
    cout << "Please input the query sequence:" << endl;
    cin >> query;

    int n = query.length();
    int m = reference.length();

    // 预处理 reference 的所有子串及其互补翻转序列
    unordered_map<string, pair<int, bool>> substring_map; // 子串 -> {起始位置, 是否翻转}
    for (int i = 0; i < m; i++) {
        string sub = "";
        for (int j = i; j < m; j++) {
            sub += reference[j];
            substring_map[sub] = {i, false}; // 标记为非互补翻转
        }
    }
    
    // 生成 reference 的互补翻转版本
    string reference_rc = get_reverse_complement(reference);
    
    // 遍历 reference_rc 的所有子串
    for (int i = 0; i < m; i++) {
        string sub = "";
        for (int j = i; j < m; j++) {
            sub += reference_rc[j];
            int len = j - i + 1;
            int original_start = m - len - i; // 转换为 reference 中的起始下标
            if (substring_map.find(sub) == substring_map.end()) {
                substring_map[sub] = {original_start, true}; // 标记为互补翻转
            }
        }
    }

    // 动态规划分割 query
    vector<int> dp(n + 1, n + 1); // 初始化为最大值
    dp[n] = 0; // 空串为 0
    // trace[i] 记录从位置 i 经过匹配后到下一个位置的信息：{next_pos, {start_pos, is_inversion}}
    vector<pair<int, pair<int, bool>>> trace(n + 1, {-1, {-1, false}});

    for (int i = n - 1; i >= 0; i--) {
        string sub = "";
        for (int j = i; j < n; j++) {
            sub += query[j];
            if (substring_map.find(sub) != substring_map.end()) {
                int next_pos = j + 1;
                if (dp[i] > 1 + dp[next_pos]) {
                    dp[i] = 1 + dp[next_pos];
                    int start = substring_map[sub].first;
                    bool is_rc = substring_map[sub].second;
                    trace[i] = {next_pos, {start, is_rc}};
                }
            }
        }
    }

    // 根据 trace 生成每个分段的结果信息，存储为 {start_pos, length, is_inversion}
    // 其中 length 为每个分段（重复单位）的长度
    vector<pair<pair<int,int>, bool>> segments;
    int pos = 0;
    while (pos < n) {
        int next = trace[pos].first;
        int start = trace[pos].second.first;
        bool is_rc = trace[pos].second.second;
        int len = next - pos;  // 重复单位的长度
        segments.push_back({{start, start + len - 1}, is_rc});
        pos = next;
    }

    // 对连续出现相同结果的 segments 进行分组合并
    // 判断条件为：起始位置和翻转标志相同。重复次数即为连续相同分段的个数，
    // 而每个分段的 length 是单位长度
    struct Output {
        int start_pos;
        int length;
        int repeat_count;
        bool is_inversion;
    };
    vector<Output> outputs;
    for (size_t i = 0; i < segments.size();) {
        int current_start = segments[i].first.first;
        // 单位长度计算为（区间右端点 - 起始 + 1）
        int unit_length = segments[i].first.second - segments[i].first.first + 1;
        bool current_rc = segments[i].second;
        int count = 1;
        size_t j = i + 1;
        while (j < segments.size() && segments[j].first.first == current_start && segments[j].second == current_rc) {
            count++;
            j++;
        }
        outputs.push_back({current_start, unit_length, count, current_rc});
        i = j;
    }

    cout << "Notice:" << endl;
    cout << "1. 本输出格式中，从前往后的顺序考察query序列中的每一段分别来自于reference的何位置，复制多少次，是否翻转等信息。"
    cout << "2. 每条记录中，start_pos 表示 query 的该段序列在 reference 里的起始位置（0 为起始索引位置），length 为重复单位长度，repeat_count 为重复次数（至少 1 次，即其本身），is_inversion 表示该段序列是否翻转。" << endl;
    
    // 输出结果，格式为：
    // {'start_pos': X, 'length': Y, 'repeat_count': Z, 'is_inversion': True/False}
    for (auto& out : outputs) {
        cout << "{'start_pos': " << out.start_pos 
             << ", 'length': " << out.length 
             << ", 'repeat_count': " << out.repeat_count 
             << ", 'is_inversion': " << (out.is_inversion ? "True" : "False") 
             << "}" << endl;
    }

    return 0;
}