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

    // 动态规划
    vector<int> dp(n + 1, n + 1); // 初始化为最大值
    dp[n] = 0; // 空串为 0
    vector<pair<int, pair<int, bool>>> trace(n + 1, {-1, {-1, false}}); // 记录转移

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

    // 输出结果
    vector<pair<pair<int, int>, bool>> result;
    int pos = 0;
    while (pos < n) {
        int next = trace[pos].first;
        int start = trace[pos].second.first;
        bool is_rc = trace[pos].second.second;
        int len = next - pos;
        result.push_back({{start, start + len - 1}, is_rc});
        pos = next;
    }

    for (auto& res : result) {
        cout << "[" << res.first.first << "," << res.first.second << "] "
             << (res.second ? "yes" : "no") << endl;
    }

    return 0;
}