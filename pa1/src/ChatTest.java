// ChatServer

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChatTest {
    private static final Charset utf8 = Charset.forName("UTF-8");

    private static final String OK = "200 OK";
    private static final String NOT_FOUND = "404 NOT FOUND";
    private static final String HTML = "text/html";
    private static final String TEXT = "text/plain";
    private static final String DEFAULT_ROOM = "DEFAULT";

    private static final Pattern PAGE_REQUEST
        = Pattern.compile("GET /([\\p{Alnum}]*/?) HTTP.*");
    private static final Pattern PULL_REQUEST
        = Pattern.compile("POST /([\\p{Alnum}]*)/?pull\\?last=([0-9]+) HTTP.*");
    private static final Pattern PUSH_REQUEST
        = Pattern.compile("POST /([\\p{Alnum}]*)/?push\\?msg=([^ ]*) HTTP.*");

    private static final String CHAT_HTML;
    static {
        try {
            CHAT_HTML = getFileAsString("../index.html");
        } catch (final IOException xx) {
            throw new Error("unable to start server", xx);
        }
    }
    
    private final String hostName = "127.0.0.1";
    private final int port;
    private final Map<String,ChatState> stateByName
        = new HashMap<String,ChatState>();
    
    class ServerThread extends Thread {
    	String roomName;
    	public void setRoomName(final String roomName){
    		this.roomName = roomName;
    	}
    	
    	public void run(){
    		try{
    			ChatTest.this.handle(roomName);
    		}
    		catch(IOException e){
    			System.err.println("Caught IOException: " + e.getMessage());
    			System.exit(-1);
    		}
    	}
    }
    
    private ServerThread[] serverThreads = new ServerThread[8];

    /**
     * Constructs a new {@link ChatServer} that will service requests
     * on the specified <code>port</code>. <code>state</code> will be
     * used to hold the current state of the chat.
     */
    public ChatTest(final int port) throws IOException {
        this.port = port;
    }

    public void runTest() throws IOException {
//        @SuppressWarnings("resource")
        for(int i = 0; i < serverThreads.length; i++){
    		String roomName = "Room"+(i/4);
    		serverThreads[i] = new ServerThread();
    		serverThreads[i].setRoomName(roomName);
    		serverThreads[i].start();
        }
    }
    
    private static String replaceEmptyWithDefaultRoom(final String room) {
    	if (room.isEmpty()) {
    		return DEFAULT_ROOM;
    	}
    	return room;
    }

    private void handle(final String roomName) throws IOException {
        try {
            System.out.println("Thread:" + Thread.currentThread());
            
            sendRequest("GET /"+roomName);
            for (int i = 1; i < 5; i++){
            	sendRequest("POST /"+roomName+"/push?msg=Thread:"+Thread.currentThread()+"Message:"+i);
            }
        } finally {

        }
    }
    
    /**
     * Writes a minimal but valid HTTP response to
     * <code>output</code>.
     */
    private void sendRequest(String content) throws IOException {
    	content += " HTTP/1.1";
    	final Socket connection = new Socket(hostName, port);
        final OutputStream xo = new BufferedOutputStream(connection.getOutputStream());
        final byte[] data = content.getBytes(utf8);
        xo.write(data);
        xo.flush();
        connection.close();

        System.out.println(Thread.currentThread() + ": send request " + content);
    }

    private synchronized ChatState getState(final String room) {
        ChatState state = stateByName.get(room);
        if (state == null) {
            state = new ChatState(room);
            stateByName.put(room, state);
        }
        return state;
    }

    /**
     * Reads the resource with the specified path as a string, and
     * then returns the string.
     */
    private static String getFileAsString(final String path)
        throws IOException {
    	byte[] fileBytes = Files.readAllBytes(Paths.get(path));
    	return new String(fileBytes, utf8);
    }

    /**
     * Runs a chat server, with a default port of 8080.
     */
    public static void main(final String[] args) throws IOException {
        final int port = args.length == 0 ? 8080 : Integer.parseInt(args[0]);
        new ChatTest(port).runTest();
    }
}
