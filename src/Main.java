import java.util.*;
public class Main {
    public static void main(String[] args){
        Scanner sc =new Scanner(System.in);
        int[] tree1 = new int[50];
        int[] tree2 =new int[50];
        int n = sc.nextInt(); // 砍了n棵树
        for(int i=0;i<n;i++){
            int temp=sc.nextInt();
            if(temp%2==0){
                tree2[temp/2-1] = 1;
            }else{
                tree1[temp/2]=1;
            }
        }
        int res1=longnum0(tree1);
        int res2=longnum0(tree2);
        int res=Math.max(res1,res2);
        System.out.println(res);
    }
    //求最长连续0
    public static int longnum0(int[] num){
        int max=0;
        int cur=0;
        for(int i=0;i<num.length;i++){
            if(num[i]==0) cur++;
            else cur=0;
            if(cur>max){
                max=cur;
            }
        }
        return max;
    }
    /*
    public static void main(String[] args){
        Scanner sc =new Scanner(System.in);
        int n=sc.nextInt();
        int[] num= new int[n];
        for(int i=0;i<n;i++){
            num[i]=sc.nextInt();
        }
        int[] dp = new int[n];
        dp[0]=1;
        for(int i=1;i<n;i++){
            int max=0;
            //遍历i之前的
            for(int j=0;j<i;j++){
                if(num[j]==num[i] || num[i]-num[j]==1){
                    if(dp[j]>max){
                        max=dp[j];
                    }
                }
                dp[i]=max+1;
            }
        }
        System.out.println(Arrays.toString(dp));
        int maxdp=0;
        for(int i=0;i<n;i++){
            if(dp[i]>maxdp){
                maxdp=dp[i];
            }
        }
        int res=n-maxdp;
        System.out.println(res);
    }*/
    /*
    public static void main(String[] args) {
       Scanner sc = new Scanner(System.in);
       String str = sc.next();
        int cnt = 1;
        for(int i=0; i<str.length()-1;i++){
            if(str.charAt(i)!=str.charAt(i+1)) cnt++;
        }
        if(cnt<str.length()-1) System.out.println(cnt+2);
        else System.out.println(str.length());

    }*/

    public static void reversestr(StringBuilder str,int left,int right){
        if(left>=right) return;
        char temp=str.charAt(left);
        str.setCharAt(left,str.charAt(right));
        str.setCharAt(right,temp);
        reversestr(str,left+1,right-1);
    }

    public static int[] maxSlidingWindow(int[] nums, int k) {
        if(nums.length == 0 || k == 0 || k>nums.length) return new int[0];
        int[] res= new int[nums.length-k+1];
        ArrayDeque<Integer> queue = new ArrayDeque<Integer>();
        for(int i=0;i<k;i++){
            while(!queue.isEmpty() && queue.getLast()<nums[i]){
                queue.removeLast();
            }
            queue.addLast(nums[i]);

        }
        res[0] = queue.getFirst();
        for(int i=k;i>nums.length;i++){
            if(queue.getFirst()==nums[i-k]){
                queue.removeFirst();
            }
            while(!queue.isEmpty() && queue.getLast()<nums[i]){
                queue.removeLast();
            }
            queue.addLast(nums[i]);
            int temp=queue.getFirst();
            System.out.println(queue.toString());
            res[i-k+1]=queue.getFirst();
        }
        return res;
    }

}

