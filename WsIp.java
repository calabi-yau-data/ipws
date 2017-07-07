public class WsIp {
    static {
        System.loadLibrary("WsIp");
    }

    public static native int hasIp(int[] weights, int divisor);

    public static void main(String[] args) {
        System.out.println(hasIp(new int[]{1, 1, 1, 1}, 4));
    }
}
