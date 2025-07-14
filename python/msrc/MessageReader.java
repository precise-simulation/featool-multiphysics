// MessageReader (Read messages with 4 byte uint prefix detemining message length.)

// J.S. Hysing 250410.
// Copyright 2013-2025 Precise Simulation Ltd.


// javac MessageReader.java
// jar cvf message_reader.jar MessageReader.class

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.InputStream;
import java.io.InterruptedIOException;
import java.io.IOException;
import java.io.StreamCorruptedException;
import java.net.Socket;
import java.net.SocketException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class MessageReader
{

    static public int nmbpre = 4;
    static public int hsize = 32;

    private Socket m_socket;
    private int m_timeout;
    private DataInputStream data_stream;
    private byte[] bsize;

    public MessageReader(Socket socket)
        throws IOException {
        m_socket = socket;
        m_timeout = m_socket.getSoTimeout();
        data_stream = new DataInputStream(m_socket.getInputStream());
        bsize = new byte[nmbpre];
    }

    public byte[] readMessage()
        throws EOFException, IOException, SocketException, StreamCorruptedException {
        return this.readMessage(0);
    }

    public byte[] readMessage(int timeout)
        throws IOException, SocketException, StreamCorruptedException {
        boolean set_timeout = timeout != m_timeout;
        try {
            if (set_timeout) { m_socket.setSoTimeout(timeout); }
            data_stream.readFully(bsize, 0, nmbpre);
            int msize = ByteBuffer.wrap(bsize).order(ByteOrder.nativeOrder()).getInt();
            if (msize < hsize) { throw new IllegalArgumentException("Invalid message size."); }

            byte[] buffer = new byte[msize - nmbpre];
            data_stream.readFully(buffer, 0, msize - nmbpre);
            return buffer;
        }

        catch (EOFException e)
            {
                byte[] ret = new byte[1]; ret[0] = -1;
                return ret;
            }

        catch (InterruptedIOException e)
            {
                // System.out.println("Timeout occured");
                return new byte[0];
            }

        finally {
            if (set_timeout) { m_socket.setSoTimeout(m_timeout); }
        }
    }
}
