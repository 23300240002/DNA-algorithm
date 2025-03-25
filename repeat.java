import java.util.*;

public class repeat {
    // 内部类用于存储子串信息
    static class SubInfo {
        int start;
        boolean inversion;
        SubInfo(int start, boolean inversion) { this.start = start; this.inversion = inversion; }
    }
    
    // 内部类用于记录动态规划trace信息
    static class Trace {
        int next;
        SubInfo info;
        Trace(int next, SubInfo info) { this.next = next; this.info = info; }
    }
    
    // 内部类用于最终输出的记录
    static class OutputRecord {
        int startPos;
        int length;
        int repeatCount;
        boolean isInversion;
        OutputRecord(int startPos, int length, int repeatCount, boolean isInversion) {
            this.startPos = startPos;
            this.length = length;
            this.repeatCount = repeatCount;
            this.isInversion = isInversion;
        }
        @Override
        public String toString() {
            return String.format("{'start_pos': %d, 'length': %d, 'repeat_count': %d, 'is_inversion': %s}",
                    startPos, length, repeatCount, (isInversion ? "True" : "False"));
        }
    }
    
    // 计算互补翻转序列
    static String compReverse(String seq) {
        StringBuilder sb = new StringBuilder();
        for (int i = seq.length() - 1; i >= 0; i--) {
            char c = seq.charAt(i);
            switch(c) {
                case 'A': sb.append('T'); break;
                case 'T': sb.append('A'); break;
                case 'C': sb.append('G'); break;
                case 'G': sb.append('C'); break;
                default: sb.append(c);
            }
        }
        return sb.toString();
    }
    
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Please input the reference sequence:");
        String refSeq = scanner.nextLine();
        System.out.println("Please input the query sequence:");
        String qrySeq = scanner.nextLine();
        
        int m = refSeq.length();
        int n = qrySeq.length();
        
        // 构建子串映射：子串 -> (起始位置, inversion标志)
        HashMap<String, SubInfo> subMap = new HashMap<>();
        for (int i = 0; i < m; i++) {
            StringBuilder sub = new StringBuilder();
            for (int j = i; j < m; j++) {
                sub.append(refSeq.charAt(j));
                subMap.put(sub.toString(), new SubInfo(i, false));
            }
        }
        
        String refRC = compReverse(refSeq);
        for (int i = 0; i < m; i++) {
            StringBuilder sub = new StringBuilder();
            for (int j = i; j < m; j++) {
                sub.append(refRC.charAt(j));
                int len = j - i + 1;
                int origStart = m - len - i;
                String key = sub.toString();
                if (!subMap.containsKey(key)) {
                    subMap.put(key, new SubInfo(origStart, true));
                }
            }
        }
        
        // 动态规划划分 qrySeq
        int[] dp = new int[n + 1];
        Arrays.fill(dp, n + 1);
        dp[n] = 0;
        Trace[] trace = new Trace[n + 1];
        
        for (int i = n - 1; i >= 0; i--) {
            StringBuilder sub = new StringBuilder();
            for (int j = i; j < n; j++) {
                sub.append(qrySeq.charAt(j));
                String key = sub.toString();
                if (subMap.containsKey(key)) {
                    if (dp[i] > 1 + dp[j + 1]) {
                        dp[i] = 1 + dp[j + 1];
                        trace[i] = new Trace(j + 1, subMap.get(key));
                    }
                }
            }
        }
        
        // 根据trace构造各段信息，每段包含：reference区间及 inversion标志
        ArrayList<Map.Entry<int[], Boolean>> segments = new ArrayList<>();
        int pos = 0;
        while (pos < n) {
            Trace t = trace[pos];
            int next = t.next;
            int segLen = next - pos;
            int segStart = t.info.start;
            segments.add(new AbstractMap.SimpleEntry<>(new int[]{segStart, segStart + segLen - 1}, t.info.inversion));
            pos = next;
        }
        
        // 合并连续相同的段（即起始位置及 inversion标志相同），统计重复次数
        ArrayList<OutputRecord> outputs = new ArrayList<>();
        for (int i = 0; i < segments.size(); ) {
            int curStart = segments.get(i).getKey()[0];
            int unitLen = segments.get(i).getKey()[1] - segments.get(i).getKey()[0] + 1;
            boolean curInv = segments.get(i).getValue();
            int count = 1;
            int j = i + 1;
            while (j < segments.size() &&
                   segments.get(j).getKey()[0] == curStart &&
                   segments.get(j).getValue() == curInv) {
                count++;
                j++;
            }
            outputs.add(new OutputRecord(curStart, unitLen, count, curInv));
            i = j;
        }
        
        System.out.println("Notice:");
        System.out.println("1. 本输出格式中，从前往后的顺序考察query序列中的每一段分别来自于reference的何位置，复制多少次，是否翻转等信息。");
        System.out.println("2. 每条记录中，start_pos 表示 query 的该段序列在 reference 里的起始位置（0 为起始索引位置），length 为重复单位长度，repeat_count 为重复次数（至少 1 次，即其本身），is_inversion 表示该段序列是否翻转。");
        
        for (OutputRecord rec : outputs) {
            System.out.println(rec);
        }
        scanner.close();
    }
}